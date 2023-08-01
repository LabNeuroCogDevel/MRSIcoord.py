import pytest
import grid


def test_read():
    "can we read in files and append empty?"
    roi_xy = grid.read_rois(["roi1", "roi2"], "test/data/roi_pos.txt")
    generated_roi = roi_xy[0]
    # in order or read: roi1, roi2, CaudL, sid3
    assert generated_roi == ["roi1", 5, 5]
    caud_l = roi_xy[2]
    print(roi_xy)
    assert caud_l == ["CaudL", 116, 124]
