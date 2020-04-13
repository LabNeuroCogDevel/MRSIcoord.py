.PHONY: test

test:
	python -m pytest test/
	# this matlab bit works the same as the other
	diff test/data/matlab/output/kspace.1.1  /Volumes/Hera/Projects/7TBrainMech/subjs/11743_20190802/slice_PFC/MRSI_roi/raw//kspace.1.1

