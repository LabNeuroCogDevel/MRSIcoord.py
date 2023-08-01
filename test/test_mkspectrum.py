import mkspectrum
def test_mkspectrum_opts(tmpdir):
    opts = {'--kspace': None,
     '--pos': True,
     '<col>': '88',
     '<row>': '112',
     '--sires': None,
     }
    opts = mkspectrum.update_args(opts)
    assert opts['sires'] == 24

def test_mkspectrum(tmpdir):
    opts = {'--kspace': None,
     '--pos': True,
     '<outdir>': tmpdir,
     '<posfile>': None,
     '<col>': '88',
     '<row>': '112',
     '<scoutres>': '216',
     '--sires': None,
     '<siarray1>': 'test/data/siarray.1.1'}
    opts = mkspectrum.update_args(opts)
    ret = mkspectrum.mkspectrum(opts)
