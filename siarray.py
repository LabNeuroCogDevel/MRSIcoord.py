import numpy as np
from numpy import real, imag, zeros
from numpy.fft import fftshift, ifft2


class SIArray:
    def __init__(self, siname: str, res=(24, 24), pts=1024, sliceno=1):
        self.fname = siname
        self.res = res
        self.rows = res[0]
        self.cols = res[1]
        self.pts = pts
        self.sliceno = sliceno
        self.data = None

        # NB - only ever tested with rows=cols
        if(self.rows != self.cols):
            raise Exception("rows!=cols is not tested")

        self.readsi()  # populate data

    def readsi(self):
        '''readsi - read siarray.1.1 file
        adapated from LoadSInD.m
        https://stackoverflow.com/questions/2146031/what-is-the-equivalent-of-fread-from-matlab-in-python
        https://stackoverflow.com/questions/44335749/read-a-float-binary-file-into-2d-arrays-in-python-and-matlab
        '''
        with open(self.fname, 'r') as fp1:
            # if not slice 1, calculate offset
            pxs = self.res[0] * self.res[1]
            offsetptr = pxs*self.pts*4*2*(self.sliceno-1)
            fp1.seek(offsetptr)
            # little-endian float32
            SI = np.fromfile(fp1, '<4f').reshape(pxs, 2*self.pts).T
        self.data = SI

    def integrateSI(self, s, e):
        return np.sum(self.data[s:e, :], 0).reshape(self.res)

    def IFFTData(self):
        kspace = zeros([self.rows, self.cols, 2*self.pts])

        # first half of dim1 is real, second half is imaginary component
        SIData = (self.data[:self.pts, :] +
                  self.data[self.pts:, :] * np.complex(0, 1)).\
            T.reshape(self.rows, self.cols, self.pts)
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


class Scout:
    def __init__(self, scout: str, res=216):
        self.fname = scout
        self.res = res
        with open(self.fname, 'r') as fp1:
            self.data = np.fromfile(fp1, '<4f').reshape(res, res).T
            # Matlab like
            # fp1 = fopen(fname,'r');
            # scout = fread(fp1,[res res],'float');

    def RegenCoor(self, pos, vo=0, ho=0, angle=0):
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
        rotm1 = [np.cos(angle), np.sin(angle),  vo]
        rotm2 = [-np.sin(angle), np.cos(angle), ho]

        rotm = np.vstack((rotm1, rotm2))
        thirdrow = np.ones(pos.shape[0])
        postemp = np.vstack((pos.T, thirdrow))
        lastrow = [0, 0, 1]
        rotmtemp = np.vstack((rotm, lastrow))
        retpos = np.linalg.lstsq(rotmtemp, postemp)
        return(retpos[0])
