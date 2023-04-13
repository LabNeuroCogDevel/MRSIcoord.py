import numpy as np
import pytest
from scipy.io import loadmat

import lcmodel
from siarray import Scout, SIArray


def belowthres(x, y, thres=10**-10):
    return abs(x - y).max() < thres


# ## test little-endian float32
# N.B. matlab's save highest precision is 16bit?
#      coarse python to compatible types
def test_read_scout():
    scout = Scout("test/data/mprage_middle.mat")
    # see generate_mat.m
    ml_s = loadmat("test/data/matlab/scout.mat")["scout"]
    assert (scout.data.astype("uint8") == ml_s).all()


def test_read_si():
    SI = SIArray("test/data/siarray.1.1")
    # see generate_mat.m
    ml_si = loadmat("test/data/matlab/si.mat")["SI"]
    assert (ml_si == SI.data.astype("double")).all()


def test_kspace():
    SI = SIArray("test/data/siarray.1.1")
    SI.IFFTData()
    # see generate_mat.m
    ml_k = loadmat("test/data/matlab/kspace.mat")["kspace"]
    assert belowthres(ml_k, SI.kspace)


def test_retpos_3():
    pos = np.array([[130, 99], [121, 94], [113, 89]])
    scout = Scout("test/data/mprage_middle.mat")
    retpos = scout.RegenCoor(pos)
    # see generate_mat.m
    ml_rp = loadmat("test/data/matlab/retpos_3.mat")["retpos"]
    # assert (retpos == ml_rp).all()
    assert belowthres(ml_rp, retpos)


def test_shiftmap():
    SI = SIArray("test/data/siarray.1.1")
    SI.IFFTData()

    shiftmat = SI.ShiftMap(vertshift=1.49464, horzshift=-1.60098)
    # see generate_mat.m
    ml_sm = loadmat("test/data/matlab/shiftmat.mat")["SHIFTMAT"]
    assert belowthres(ml_sm, shiftmat)


def test_littleindian(tmpdir):
    SI = SIArray("test/data/siarray.1.1")
    SI.IFFTData()
    SI.savekspace(tmpdir.join("temp_kspace.1.1"), reload=True)
    with open("test/data/matlab/kspace.1.1", "r") as f:
        mlk = np.fromfile(f, "<4f").reshape(SI.rows, SI.cols, SI.pts * 2)
    assert belowthres(SI.kspace, mlk, 10**-2)
    assert belowthres(SI.kspace, mlk, 10**-6)


def test_spatialtransform():
    SI = SIArray("test/data/siarray.1.1")
    st = SI.SpatialTransform2D(vertshift=1.49464, horzshift=-1.60098)
    # see generate_mat.m
    ml_st = loadmat("test/data/matlab/spatialtransform2d.mat")["SI"]
    assert belowthres(ml_st, st, 10**-2)  # PROBLEM?!

    # use same input as matlab
    with open("test/data/matlab/kspace.1.1", "r") as f:
        kspace = np.fromfile(f, "<4f")
    SI.kspace = kspace.reshape(24, 24, 1024 * 2)
    readk_st = SI.SpatialTransform2D(vertshift=1.49464, horzshift=-1.60098)
    assert belowthres(readk_st, ml_st)  # FIX?


def test_posshift():
    scout = Scout("test/data/mprage_middle.mat")
    SI = SIArray("test/data/siarray.1.1")
    pos = np.array([[130, 99], [121, 94], [113, 89]])
    pospp = SI.pos_shift(scout, pos)
    ml_pp = loadmat("test/data/matlab/pospp.mat")["poslpp"]
    assert belowthres(pospp, ml_pp)


def test_spectrum():
    pos = np.array([[130, 99], [121, 94], [113, 89]])
    scout = Scout("test/data/mprage_middle.mat")
    SI = SIArray("test/data/siarray.1.1")
    # use the same kspace data
    with open("test/data/matlab/kspace.1.1", "r") as f:
        SI.kspace = np.fromfile(f, "<4f")

    spectrums = SI.ReconCoordinates3(scout, pos)
    ml_s3 = loadmat("test/data/matlab/spectrum_113.89")["spectrum"]
    assert belowthres(spectrums[2, :], ml_s3)


def test_spectrum_save(tmpdir):
    pos = np.array([[130, 99], [121, 94], [113, 89]])
    scout = Scout("test/data/mprage_middle.mat")
    SI = SIArray("test/data/siarray.1.1")
    SI.IFFTData()
    SI.savekspace(tmpdir.join("temp_kspace.1.1"), reload=True)

    spectrums = SI.ReconCoordinates3(scout, pos)
    ml_s3 = loadmat("test/data/matlab/spectrum_113.89")["spectrum"]
    assert belowthres(spectrums[2, :], ml_s3, 10**-6)


def test_1d_ifft():
    ml = np.loadtxt("test/data/csi.raw.112.88", skiprows=7)
    py = lcmodel.ifft_spec("test/data/spectrum.112.88", npoints=1024)
    py2 = np.stack([py.real, py.imag]).T
    np.testing.assert_allclose(ml, py2, atol=6e-5)


def test_1d_ifft_csiraw_header():
    lcm = lcmodel.LCModel("test/data/spectrum.112.88", lcmodel_path=None)
    pycsiraw = lcm.write_raw_jref(csi_raw_fname=None)
    pycsiraw = pycsiraw.replace("ID='None'", "ID='csi.raw'")
    with open("test/data/csi.raw.112.88", "r") as f:
        mlcsiraw = f.read()

    # only check the header: first 7 lines
    n = 8
    assert pycsiraw.split("\n")[0:n] == mlcsiraw.split("\n")[0:n]

    # python -m pytest test/test_match_matlab.py -k 1difft -vv
    # str diff on very large strings is slow!

    # E         -    2.658653E-01  -1.907578E-01
    # E         ?          --             --
    # E         +    2.658546E-01  -1.907615E-01
    # E         ?         ++             ++
    # E         -   -2.113438E-01  -3.969765E-02
    # E         ?         ^ -           ^^^^
    # E         +   -2.113213E-01  -3.971634E-02
    # E         ?         ^^           ++ ^^
    # E         -    1.346388E-01  -6.609750E-02
    # E         +    1.346450E-01  -6.610133E-02
    # E         -   -3.491483E-01  -8.017540E-03
    # E         ?          ^^       ^ ^^  --
    # E         +   -3.491475E-01  -7.996765E-03
    # E         ?          ^^       ^ ^^^ +
