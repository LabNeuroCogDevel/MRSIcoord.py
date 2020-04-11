import numpy as np
from scipy.io import loadmat
from siarray import Scout, SIArray
import pytest

def belowthres(x,y,thres=10**-10):
    return(abs(x - y).max() < thres)

# ## test little-endian float32
# N.B. matlab's save highest precision is 16bit?
#      coarse python to compatible types
def test_read_scout():
    scout = Scout('test/data/mprage_middle.mat')
    # see generate_mat.m
    ml_s = loadmat('test/data/matlab/scout.mat')['scout']
    assert (scout.data.astype('uint8') == ml_s).all()


def test_read_si():
    SI = SIArray('test/data/siarray.1.1')
    # see generate_mat.m
    ml_si = loadmat('test/data/matlab/si.mat')['SI']
    assert (ml_si == SI.data.astype('double')).all()


def test_kspace():
    SI = SIArray('test/data/siarray.1.1')
    SI.IFFTData()
    # see generate_mat.m
    ml_k = loadmat('test/data/matlab/kspace.mat')['kspace']
    assert belowthres(ml_k, SI.kspace)

def test_retpos_3():
    pos = np.array([[130, 99], [121, 94], [113, 89]])
    scout = Scout('test/data/mprage_middle.mat')
    retpos = scout.RegenCoor(pos)
    # see generate_mat.m
    ml_rp = loadmat('test/data/matlab/retpos_3.mat')['retpos']
    #assert (retpos == ml_rp).all()
    assert belowthres(ml_rp, retpos)

def test_shiftmap():
    SI = SIArray('test/data/siarray.1.1')
    SI.IFFTData()
    
    shiftmat = SI.ShiftMap(vertshift=1.49464, horzshift=-1.60098)
    # see generate_mat.m
    ml_sm = loadmat('test/data/matlab/shiftmat.mat')['SHIFTMAT']
    assert belowthres(ml_sm, shiftmat)

def test_spatialtransform():
    SI = SIArray('test/data/siarray.1.1')
    st = SI.SpatialTransform2D(vertshift=1.49464, horzshift=-1.60098)
    # see generate_mat.m
    ml_st = loadmat('test/data/matlab/spatialtransform2d.mat')['SI']
    assert belowthres(ml_st, st, 10**-2) # PROBLEM?!
    
    # use same input as matlab
    with open('test/data/matlab/kspace.1.1', 'r') as f:
        kspace = np.fromfile(f,'<4f')
    SI.kspace = kspace.reshape(24,24,1024*2)
    readk_st = SI.SpatialTransform2D(vertshift=1.49464, horzshift=-1.60098)
    assert belowthres(readk_st, ml_st) # FIX?

def test_posshift():
    scout = Scout('test/data/mprage_middle.mat')
    SI = SIArray('test/data/siarray.1.1')
    pos = np.array([[130,99], [121, 94], [113, 89]])
    pospp = SI.pos_shift(scout, pos)
    ml_pp = loadmat('test/data/matlab/pospp.mat')['poslpp']
    assert belowthres(pospp, ml_pp)

def test_spectrum():
    pos = np.array([[130,99], [121, 94], [113, 89]])
    scout = Scout('test/data/mprage_middle.mat')
    SI = SIArray('test/data/siarray.1.1')
    # use the same kspace data
    with open('test/data/matlab/kspace.1.1', 'r') as f:
        SI.kspace = np.fromfile(f, '<4f')

    spectrums = SI.ReconCoordinates3(scout, pos)
    ml_s3 = loadmat('test/data/matlab/spectrum_113.89')['spectrum']
    assert belowthres(spectrums[2,:], ml_s3)
    
