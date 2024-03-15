#!/usr/bin/env python
"""
wrap around lcmodel.
spectrum data extracted from siarray.1.1 (1024 complex points): IFFT around coordinate
see mkspectrum (siarray.py)
"""
import os
import subprocess

import numpy as np
import numpy.typing as npt
from numpy.fft import ifft

from siarray import SIArray


def write_string(fname: str, contents: str) -> None:
    """
    write string to file. do nothing if fname is None or empty.
    (empty fname used for testing)
    """
    if not fname:
        return
    with open(fname, "w") as f:
        f.write(contents)


def ifft_spec(spec_fname: str, npoints: int) -> npt.NDArray[complex]:
    """
    read in the and ifft complex spectrum of at individual coordinate
    output used for LCModel.write_raw_jref

    :param spec_fname: little endian float input vector file
    :param npoints: number of complex elements in input (1/2 total: real + imag)
    :returns: 1D ifft complex values
    """

    si = SIArray(spec_fname, res=(1, 1), pts=npoints)
    SP_SLS = si.to_complex().squeeze()
    TD_SLS = ifft(SP_SLS)
    # 2nd and every other after needs to have sign flipped
    TD_SLS[1::2] = -1 * TD_SLS[1::2]
    return TD_SLS


class LCModel:
    def __init__(
        self,
        spec_file: str,
        lcmodel_path: str,
        npoints=1024,
        axPPM=297.211197,
        RBW=3000.0,
        TE=17,
        basis="lcmodel/basis-sets/gamma_7TJref_te34_297mhz_no2HG_1.basis",
    ):
        """
        LCModel object from a spectrum file (list of complex values from placed coordinate)
        :param spec_file: file w/ single vector of complex values.
                          little endian. see SIArray for reading notes
        :param lcmodel_path: path to lcmodel binary
        :param npoints:   number of complex values expected in spec file
        :param axPPM:     7T default, prisma value=123.254000
        :param RBW:       Hz, 7T default 3T value=1301
        :param TE:        Echo Time. 7T default (17ms)
        :param basis:     basis function file for lcmodel control input
        :returns: LCModel object
        """
        self.axPPM = axPPM
        self.RBW = RBW
        self.TE = TE
        self.basis = os.path.abspath(basis)
        self.npoints = npoints
        if not lcmodel_path:
            lcmodel_path = os.path.join(os.path.dirname(__file__), "lcmodel/lcmodel")
        self.lcmodel_path = os.path.abspath(lcmodel_path)
        self.spec_file = spec_file  # only used to create default directory
        self.complex_ifft = ifft_spec(spec_file, self.npoints)

        # updated by write_raw_jref, used by write_control
        self.csi_raw_fname = None
        # updated by write_control used by lcmodel()
        self.control_fname = None

    def run(self, outdir=None, redo=False):
        """
        create lcmodel temporary files and run
        """
        if not outdir:
            outdir = self.spec_file + ".dir"

        output_sheet = os.path.join(outdir,'spreadsheet.csv')
        if not redo and os.path.isfile(output_sheet):
            print(f"already have '{output_sheet}'. use lcmodel.run(...,redo=True) in python to rerun")
            return 0

        os.makedirs(outdir, exist_ok=True)
        pwd = os.getcwd()
        os.chdir(outdir)
        raw_str = self.write_raw_jref()
        control_str = self.write_control()
        ret = self.lcmodel()
        os.chdir(pwd)
        return ret

    def lcmodel(self):
        """run lcmodel on control file. be sure to create temporary files first."""
        if not os.access(self.lcmodel_path, os.X_OK) or not os.path.isfile(self.lcmodel_path):
            raise Exception(f"do not have valid lcmodel binary '{self.lcmodel_path}'")

        cmd = f"{self.lcmodel_path} < {os.getcwd()}/{self.control_fname}"
        return subprocess.run(cmd, shell=True).returncode
        # print(cmd)
        # return os.system(cmd)

    def write_raw_jref(self, csi_raw_fname="csi.raw") -> str:
        """
        csi.raw lcmodel input file. adapted from MRRC Matlab code.
        simplified to always use J-REF
        NB. 15.6E formatting specific to lcmodel on linux?

        :param fname: file to write to
                      if None, will only return file contents as string
        :returns: string of file csi.raw contents
        """

        # header is mostly static text (dynamic for fname, te, and ax).
        # spacing, additional . after TE matches ml output. unclear if needed
        raw = f""" $SEQPAR
 ECHOT ={2*self.TE}.
 HZPPPM={self.axPPM}
 SEQ = 'J-REF'
 $END
 $NMID ID='{csi_raw_fname}', FMTDAT='(2e15.6)'
 TRAMP= 1.0, VOLUME=1 $END\n"""
        # 15 characters with 6 decimal places and Exx
        raw += "\n".join(
            [f"{np.real(x):15.6E}{np.imag(x):15.6E}" for x in self.complex_ifft]
        )
        raw += "\n"  # to match matlab, not necessarily b/c lcmodel requires

        # can test function with fName=None to just check string
        write_string(csi_raw_fname, raw)
        self.csi_raw_fname = csi_raw_fname
        return raw

    def write_control(self, control_fname="csi.control") -> str:
        """
        construct lcmodel control file
        """

        self.control_fname = control_fname

        # .0003333 (no leading zero)
        DELTAT = str(round(1 / self.RBW, 7)).replace("0.", ".")

        control_str = f""" $LCMODL
 TITLE= '{os.getcwd()}/{control_fname}' 
 FILBAS= '{self.basis}'
 FILRAW='{self.csi_raw_fname}' 
 FILCOO='csi.coord' 
 FILCSV='spreadsheet.csv' 
 FILPS='csi.ps' 
 HZPPPM={self.axPPM}
 NUNFIL={self.npoints}
 DELTAT={DELTAT}
 NDCOLS = 1
 NDROWS = 1
 NDSLIC = 1
 ICOLST = 1
 ICOLEN = 1
 IROWST = 1
 IROWEN = 1
 ISLICE = 1
 neach = 10
 nameac(10) = 'mI' 
 nameac(9) = 'Gln' 
 nameac(8) = 'GABA' 
 nameac(7) = 'GSH' 
 nameac(6) = 'Glu' 
 nameac(5) = 'GPC' 
 nameac(4) = 'NAAG' 
 nameac(3) = 'NAA' 
 nameac(2) = 'Cho' 
 nameac(1) = 'Cre' 
 IPAGE2 = 0
 PPMST=4
 PPMEND=1.8
 LCOORD=9
 LCSV=11
 $END"""
        write_string(control_fname, control_str)
        return control_str


def run_lcmodel(spectrum_files):
    lcmodel_path = os.path.join(os.path.dirname(__file__), "lcmodel/lcmodel")
    basis_path = os.path.join(
        os.path.dirname(lcmodel_path),
        "basis-sets/gamma_7TJref_te34_297mhz_no2HG_1.basis",
    )
    for spec_file in spectrum_files:  # "test/data/spectrum.112.88"
        print(f"running lcmodel for: '{spec_file}'")
        lcm = LCModel(spec_file, lcmodel_path, basis=basis_path)
        lcm.run()


if __name__ == "__main__":
    import sys

    run_lcmodel(sys.argv[1:])
