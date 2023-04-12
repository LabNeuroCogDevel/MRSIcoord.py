"""
wrap around lcmodel.
spectrum data extracted from siarray.1.1 (1024 complex points): IFFT around coordinate
see mkspectrum (siarray.py)
"""
from siarray import SIArray
import numpy as np
from numpy.fft import ifft

def ifft_spec_WRONG(spec_fname):
    """spectrum.xx.yy is same format as siarray.1.1. we can reuse fft from there (TODO: confirm)
    @param spec_fname - spectrum.xx.yy file output from coordinate placement
    @return kspace (IFFT) 1024 complex values"""
    si = SIArray(spec_fname, res=(1, 1))
    si.IFFTData()
    arr = si.kspace.squeeze()
    # pair complex data
    n = len(arr)//2 # 1024
    #paired = np.stack([arr[0:n], arr[n:]])
    paired = np.stack([arr[0::2], arr[1::2]])
    return paired

def ifft_spec(spec_fname, Nkeep=1024, Nx=1,Ny=1,NslsProc=1):
    """
    read in the and ifft the (likely) 1x1024 complex specturm saved by coordinate placement
    output used for write_raw_jref
    """

    # si = SIArray('out/spectrum.112.88', res=(1, 1))
    # SL = si.data
    # SP_SLS = si.to_complex()

    SL = np.fromfile(spec_fname, "<f").reshape(2 * Nkeep,Nx,Ny)
    # into complex
    SP_SLS = np.zeros([Nkeep,Nx,Ny,NslsProc], dtype='complex') # f is float32 => single-precision
    real = SL[0:Nkeep,:,:]
    imag = 1j * SL[Nkeep+0:2*Nkeep,:,:]
    SP_SLS[:,:,:,0] = real + imag
    
    SP_SLS=np.squeeze(SP_SLS)
    TD_SLS = ifft(SP_SLS)
    TD_SLS[1::2] = (-1) * TD_SLS[1::2]
    return TD_SLS
    

def write_raw_jref(fName, complex_ifft, TE=17, axPPM=297.211197):
    """
    csi.raw lcmodel input file. adapated from MRRC Matlab code. simplified to always use J-REF
    @fName file to write to
    @dataLC data from spectrum IFFT at corrdinate. expect 1024 complex values for 7T
    @TE echo time (default to 17ms)
    @axPPM ???? (default to 7T value, matlab mentions Prisma value=123254000)
    
    NB. 15.6E formating specific to lcmodel on linux?
    """
    # header is mostly static text (dynamic for fname, te, and ax).
    # spacing, additional . after TE matches ml output. unclear if needed
    raw=f""" $SEQPAR
 ECHOT ={2*TE}.
 HZPPPM={axPPM}
 SEQ = 'J-REF'
 $END
 $NMID ID='{fName}', FMTDAT='(2e15.6)'
 TRAMP= 1.0, VOLUME=1 $END\n"""
    # 15 characters with 6 decimal places and Exx
    raw+="\n".join([f"{np.real(x):15.6E}{np.imag(x):15.6E}" for x in complex_ifft])
    raw+="\n" # to match matlab, not necessarily b/c lcmodel requires

    # can test function with fName=None to just check string
    if fName:
        with open(fName,'w') as f:
            print(f, raw)
    return raw

