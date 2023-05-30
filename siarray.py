"""
Port of MRRC's sid3 and SVR1HFinal (matlab)

rhea:/opt/ni_tools/matlab_toolboxes/MRRC/
"""
import numpy as np
import numpy.typing as npt
from numpy import imag, real, zeros
from numpy.fft import fftshift, ifft2


# see DisplaySIImage.m,  AdjustSIImage.m, ReadSliderText2.m
def adjust(img, brightness, contrast, gamma=1):
    """[a,b] to [c,d] with gamma. like imadjust in matlab"""
    a = np.min(img)
    b = np.min(img)
    c = brightness
    d = (1.0-c)*contrast
    # imadjust
    return (((img - a) / (b - a)) ** gamma) * (d - c) + c


class Offsets:
    def __init__(self, vo=0, ho=0, angle=0):
        """offsets
        @param vo vertical offset
        @param ho horizons offset
        @param angle in radians
        """
        self.vo = vo
        self.ho = ho
        self.angle = angle
        if vo != 0 or ho != 0 or angle != 0:
            raise Warning("untested offsets/angle rotation!")
        rotm1 = [np.cos(angle), np.sin(angle), vo]
        rotm2 = [-np.sin(angle), np.cos(angle), ho]
        self.rotm = np.vstack((rotm1, rotm2))


class Scout:
    def __init__(self, scout: str, res=216):
        self.res = res
        self.fname = scout
        self.data = None
        if scout:
            with open(self.fname, "r") as fp1:
                self.data = np.fromfile(fp1, "<4f").reshape(res, res).T
                # Matlab like
                # fp1 = fopen(fname,'r');
                # scout = fread(fp1,[res res],'float');

    def RegenCoor(self, pos, offsets=Offsets()):
        """
        optional rotate coords (default to not) and put position into scout space
        @param pos like [[row, column],...]
        @return retpos
        >>> pos = np.array([[130,99], [121, 94]])
        """

        # Note from MATLAB code:
        #  `pos` are in SID orientation, but in Transfmd/cropped image-space
        #  purely for user convenience
        #  so since reconstruction is still done in MATLAB, need MATLAB orientation
        res = self.res
        pos = res + 2 - pos

        # create transformation matrix, 2rows x 3columns
        rotm = offsets.rotm

        thirdrow = np.ones(pos.shape[0])
        postemp = np.vstack((pos.T, thirdrow))
        lastrow = [0, 0, 1]
        rotmtemp = np.vstack((rotm, lastrow))
        retpos = np.linalg.lstsq(rotmtemp, postemp, rcond=None)
        return np.array(retpos[0])


class Shifts:
    """class to easily pass around default shift settings"""

    def __init__(self, shiftvolume=1, rotangle=0):
        # toggleon=1, hanningon=1,
        # flipvert=0, fliphorz=0, flipslices=0,
        self.shiftvolume = shiftvolume
        self.rotangle = rotangle


class SIArray:
    def __init__(self, siname: str, res=(24, 24), pts=1024, sliceno=1, shift=Shifts()):
        self.fname = siname
        self.res = res
        self.rows = res[0]
        self.cols = res[1]
        self.pts = pts
        self.sliceno = sliceno
        self.data = None
        self.kspace = None
        self.shift = shift

        # NB - only ever tested with rows=cols
        if self.rows != self.cols:
            raise Exception("rows!=cols is not tested")

        self.readsi()  # populate data

    def readsi(self):
        """readsi - read siarray.1.1 file
        adapated from LoadSInD.m
        https://stackoverflow.com/questions/2146031/what-is-the-equivalent-of-fread-from-matlab-in-python
        https://stackoverflow.com/questions/44335749/read-a-float-binary-file-into-2d-arrays-in-python-and-matlab
        """
        with open(self.fname, "r") as fp1:
            # if not slice 1, calculate offset
            pxs = self.res[0] * self.res[1]
            offsetptr = pxs * self.pts * 4 * 2 * (self.sliceno - 1)
            fp1.seek(offsetptr)
            # little-endian float32
            SI = np.fromfile(fp1, "<4f").reshape(pxs, 2 * self.pts).T
        self.data = SI

    def integrateSI(self, s, e=None):
        """sum from start to end. sid3:IntegrateSI.m"""

        if not e:
            e = self.data.shape[0]
        return np.sum(self.data[s:e, :], 0).reshape(self.res)

    def to_complex(self) -> npt.NDArray[complex]:
        """
        first half of self.data is real. second half is imag
        combine as complex
        """
        real = self.data[: self.pts, :]
        imag = self.data[self.pts :, :] * complex(0, 1)
        return (real + imag).T.reshape(self.rows, self.cols, self.pts)

    def IFFTData(self):
        """inverse FFT to get back to kspace"""
        kspace = zeros([self.rows, self.cols, 2 * self.pts])

        # first half of dim1 is real, second half is imaginary component
        SIData = self.to_complex()
        # from matlab:
        #  SIData(15,13,82) == -1.0690e+02 - 5.3431e+01i
        #  SIData(20,6,507) == 2.2938e+02 + 1.8219e+02i
        # matches
        # SIData[14,12,81]
        # SIData[19,5,506]

        # ifft spatial a la MRRC matlab code
        for a in range(self.pts):
            temp = fftshift(ifft2(SIData[:, :, a]))
            kspace[:, :, a] = real(temp)
            kspace[:, :, self.pts + a] = imag(temp)

        self.kspace = kspace

    def savekspace(self, fname, reload=False):
        """save kspace to little indian encoded file
        optionally, reload it for lossy matlab match

        also note [ 1 2; 3 4]
          Matlab=column-major (1 3 2 4),
          python=row-major    (1 2 3 4)
        https://scottstaniewicz.com/articles/python-matlab-binary/
        >>> SI.kspace[0,0,0].tobytes()               == b'j\x1c\xf7oVk\x98?'
        >>> SI.kspace.astype('<f4')[0,0,0].tobytes() == b'\xb3Z\xc3<'
        """
        k = self.kspace.astype("<f4")
        with open(fname, "wb") as f:
            k.tofile(f)
        with open(fname, "r") as f:
            # little-endian float
            kn = np.fromfile(f, "<4f").reshape(self.rows, self.cols, self.pts * 2)
        if reload:
            self.kspace = kn
        return kn

    def ShiftMap(self, vertshift, horzshift):
        """
        @param shift - how to manipulate
        @return SHIFTMAT (N.B. transpose at the end to match matlab)
        """
        if not self.shift.shiftvolume:
            raise Warning("UNTESTED no shift")
            SHIFTMAT = np.ones((self.rows, self.cols)) + complex(0, 0)
            return SHIFTMAT.T

        # as saved by kspace.1.1
        r = (np.arange(self.rows) - self.rows / 2) * horzshift / self.rows
        c = (np.arange(self.cols) - self.cols / 2) * vertshift / self.cols
        rr, cc = np.meshgrid(r, c)
        angle = (rr + cc) * 2 * np.pi
        SHIFTMAT = np.exp(angle * complex(0, 1))
        return SHIFTMAT.T

    def pos_shift(self, scout: Scout, pos: np.ndarray):
        """
        get shift valuse from t1 in scout space row col pair
        @param scout - scout used for res, pasted to scout.RegenCorr
        @param pos - [ [row,col], [row,col], ... ]
        @return pospp - shifted positions
        """
        # row, col
        # in ml code there is posl and posr
        # postemp has a new bottom row of all ones
        postemp = scout.RegenCoor(pos)

        # input('These below are the SID 256x256 orientation coords!');
        posadj = scout.res + 2 - postemp

        # convert the offsets to fractional pixel shifts
        # transpose, drop 3rd column (was 3rd row).
        #  adjust by csi.res/scout.res
        fracpx = (
            (posadj.T[:, :2] - 1 - scout.res / 2)
            * np.array([self.rows, self.cols])
            / scout.res
        )
        # half px shift in both r and c dxns
        pospp = -1 * fracpx + 0.5
        return pospp

    def SpatialTransform2D(self, vertshift, horzshift):
        """
        @param vertshift
        @param horzshift
        @param shift  rotangle and other settings
        @return transformed kspace about position
        """
        # data look like matlab
        # kspSI.shape == (2048, 576)
        if self.kspace is None:
            self.IFFTData()
        SHIFTMAP = self.ShiftMap(vertshift, horzshift)

        kspSI = self.kspace.reshape(self.rows * self.cols, self.pts * 2).T
        # convert to complex - first half is read, second is half is complex
        data1 = kspSI[: self.pts, :] + kspSI[self.pts :, :] * complex(0, 1)
        # make 3d
        data_3d = data1.reshape(self.pts, self.rows, self.cols).T
        shifted_3d = np.zeros(data_3d.shape, dtype="complex128")

        for a in np.arange(self.pts):
            # voxel shift the data
            data2 = data_3d[:, :, a] * SHIFTMAP
            data2 = np.fft.fft2(data2)

            # rotate images - code from matlab. not implemented (tested)
            if self.shift.rotangle:
                raise Exception("UNTESTED/UNIMPLEMENTED rotation")
                # data2 = imrotate(data2,rotangle,'nearest','crop');
                # maybe?
                # data2c = scipy.misc.imrotate(data2,shift.rotangle,'nearest')
            shifted_3d[:, :, a] = data2

        # put back in 2d like ml code. second data1 in matlab code
        # that is. stack the real and imaginary
        kspSI_shift = shifted_3d.T.reshape(self.pts, self.rows * self.cols)
        ret = np.vstack((np.real(kspSI_shift), np.imag(kspSI_shift)))
        return ret

    def spectrum(self, vertshift, horzshift):
        """
        @param vertshift, horzshift - row, col after pos_shift
        @return spectrum pts*2 array real then complex at shifted position
        """
        st = self.SpatialTransform2D(vertshift, horzshift)
        # 1 row past center? # 300 if 24x24
        extract_pos = int((self.rows + 1) * self.cols / 2)

        # get first and second half of 2048 length vector (for 24x24 siarray)
        srd = st[: self.pts, extract_pos]
        sid = st[self.pts :, extract_pos]
        # commented out in ML code too
        # totalshift = (vertshift+horzshift)*3.14159
        totalshift = 0
        scd = (srd + sid * complex(0, 1)) * np.exp(complex(0, 1) * totalshift)
        spectrum = np.concatenate((np.real(scd), np.imag(scd)))
        return spectrum

    def ReconCoordinates3(
        self, scout: Scout, pos, writedir=None, specprefix="spectrum"
    ):
        """
        generate spectrum from a given rorig coordinate
        @param pos - row, col postions to generate spectrum
        @return (spectrums,files) - n_pos x self.pts*2 array of spectrums and output filenames
        >>> pos = np.array([[130, 99], [121, 94], [113, 89]])
        """

        # row, col
        # in ml code there is posl and posr
        numrecon = pos.shape[0]
        pospp = self.pos_shift(scout, pos)

        # load the fractional pixel shifts
        # do the reconstruction
        # convert the data matrix to appropriate formt
        spectrums = zeros((numrecon, self.pts * 2), dtype="float64")
        filenames = [None] * numrecon
        for m in np.arange(numrecon):
            spectrum = self.spectrum(pospp[m, 0], pospp[m, 1])
            spectrums[m, :] = spectrum

            if writedir:
                row = pos[m, 0]
                col = pos[m, 1]
                filenames[m] = "%s/%s.%d.%d" % (writedir, specprefix, row, col)
                with open(filenames[m], "wb") as f:
                    spectrum.astype("<f4").tofile(f)

        return (spectrums,filenames)


def ignored_regencor_scoutarray2(scout: Scout, rotm):
    """this was part of RegenCoor, but is not used?"""
    res = scout.res
    # ## CROPPING VALUES
    # generate the rotated and translated scout again
    # this whole section to calculate cropping values
    scoutarray = scout.data.T.reshape(1, res * res).flatten()
    # point-by-point matrix of loci, 3 rows x res*res columns=basm=scoutmatrx
    seq = np.arange(res) + 1
    basm = np.vstack((np.tile(seq, res), np.repeat(seq, res), np.ones(res**2)))
    scoutmatrx = basm
    # rotated imagearray (2r x 3c) * (3r x 65536 cols)
    xx = np.dot(rotm, scoutmatrx)
    # 2rows by 65536cols matrix ea column (row locus,col locus)
    xx = round(xx)
    scoutarray2 = zeros((res, res))
    # STILL ML CODE
    # assign the intensities accdg to scout and scoutarray
    # for aa=1:res*res
    #     if ((xx(1,aa)>=1) & (xx(2,aa)>=1))
    #         if ((xx(1,aa)<=256) & (xx(2,aa)<=256))
    #         scoutarray2(xx(1,aa),xx(2,aa)) = scoutarray(aa)
