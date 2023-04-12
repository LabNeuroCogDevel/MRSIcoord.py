import pytest
import grid

def test_read():
    roi_xy = grid.read_rois(["roi4","roi5"], "test/data/roi_pos.txt")
    assert(roi_xy[0] == ['roi4', 5, 5])
    assert(roi_xy[2] == ['roi1', 112.0, 88.0])
