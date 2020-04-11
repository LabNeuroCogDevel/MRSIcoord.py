import numpy as np
from scipy.io import loadmat
from siarray import Scout, SIArray

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


def test_retpos():
    pos = np.array([[130, 99], [121, 94]])
    scout = Scout('test/data/mprage_middle.mat')
    retpos = scout.RegenCoor(pos)
    # see generate_mat.m
    ml_rp = loadmat('test/data/matlab/retpos.mat')['retpos']
    assert (retpos == ml_rp).all()

def test_shiftmap():
    SI = SIArray('test/data/siarray.1.1')
    SI.IFFTData()
    shiftmat = SI.ShiftMap()
    # see generate_mat.m
    ml_sm = loadmat('test/data/matlab/shiftmat.mat')['SHIFTMAT']
    assert belowthres(ml_sm, shiftmat)

def test_spatialtransform():
    pass
    SI = SIArray('test/data/siarray.1.1')
    st = SI.SpatialTransform2D()
    # see generate_mat.m
    ml_st = loadmat('test/data/matlab/spatialtransform2d.mat')['SI']
    assert belowthres(ml_st, st, 10**-2) # PROBLEM?!
