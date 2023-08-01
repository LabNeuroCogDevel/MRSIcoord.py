.PHONY: test gui_example

test:
	python -m pytest test/
	# this matlab bit works the same as the other
	diff test/data/matlab/output/kspace.1.1  /Volumes/Hera/Projects/7TBrainMech/subjs/11743_20190802/slice_PFC/MRSI_roi/raw//kspace.1.1

## compile lcmodel from source (via github) and copy here
lcmodel/lcmodel:
	git clone https://github.com/schorschinho/LCModel lcmodel-fortran
	make -C lcmodel-fortran
	cp lcmodel-fortran/binaries/linux/lcmodel $@

## GUI program
gui_example:
	./grid.py  -s test/data/siarray.1.1  -r test/data/rorig.nii  -i test/data/roi_pos.txt  --rois roi4 roi6

## Standalone programs

# example mkspectrum call
# also see test/data/pos.txt
# x,y 116,124 => 216+1- => 101,93 => sid3: 93,101
# L Caud in /Volumes/Hera/Projects/7TBrainMech/subjs//11743_20190802/slice_PFC/MRSI_roi/13MP20200207/LT
# diff -y <(datamash -t, transpose < out/spectrum.93.101.dir/spreadsheet.csv ) <(datamash -t, transpose < /Volumes/Hera/Projects/7TBrainMech/subjs//11743_20190802/slice_PFC/MRSI_roi/LCModel/v2idxfix/spectrum.93.101.dir/spreadsheet.csv )

out/spectrum.93.101.py: test/data/siarray.1.1
	@mkdir -p out
	./mkspectrum test/data/siarray.1.1 216 --pos 93 101 out/
	mv out/spectrum.93.101 $@

out/spectrum.93.101.py.dir/: out/spectrum.93.101.py
	./lcmodel.py $<

## MATLAB (octave) vesion
out/spectrum.93.101.ml: test/data/siarray.1.1
	octave --eval "addpath('$$PWD/matlab'); gen_spectrum('test/data/siarray.1.1', 216, [93 101],'out/')"
	mv out/spectrum.93.101 $@

out/spectrum.93.101.ml.dir/: out/spectrum.93.101.ml
	octave --eval "addpath('$$PWD/matlab/lcmodel/'); cd out/; lcmodel_spec('spectrum.93.101.ml')"
