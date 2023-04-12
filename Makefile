.PHONY: test gui_example

test:
	python -m pytest test/
	# this matlab bit works the same as the other
	diff test/data/matlab/output/kspace.1.1  /Volumes/Hera/Projects/7TBrainMech/subjs/11743_20190802/slice_PFC/MRSI_roi/raw//kspace.1.1

# example mkspectrum call
# also see test/data/pos.txt
out/spectrum.112.88:
	@mkdir -p out
	./mkspectrum test/data/siarray.1.1 216 --pos 112 88 out/

gui_example:
	./grid.py  -s test/data/siarray.1.1  -r test/data/rorig.nii  -i test/data/roi_pos.txt  --rois roi4 roi6
