import lcmodel
import pytest
import os.path

def test_lcmodel(tmpdir):
    lcmodel_path="lcmodel/lcmodel"
    if not os.path.isfile(lcmodel_path):
        pytest.skip("no lcmodel binary")

    lcm = lcmodel.LCModel("test/data/spectrum.112.88", lcmodel_path)
    ret = lcm.run(tmpdir)
    assert os.path.isfile(os.path.join(tmpdir,"spreadsheet.csv"))
    assert ret == 0
