from scipy.io import loadmat
from siarray import Scout, SIArray


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
    assert abs(ml_k - SI.kspace).max() < 10**-10
