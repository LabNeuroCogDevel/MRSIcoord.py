.PHONY: test gui_example depends

test: out/test_py.txt
test_all: out/test_results.txt

.venv:
	virtualenv .venv/
	source .venv/bin/activate && \
	pip install -r requirements.txt

out/test_results.txt: out/test_ml.txt out/test_py.txt
	echo cat $^ > $@

out/test_ml.txt: test/data/11743_20190802_kspace.1.1 test/data/matlab/kspace.1.1 lcmodel/lcmodel | out
	# check this matlab bit works the same as the original MRRC version
	diff test/data/11743_20190802_kspace.1.1 test/data/kspace.1.1  | tee -a $@

out/test_py.txt: $(wildcard test/*py) lcmodel/lcmodel test/data/matlab/kspace.1.1  | out
	# python code does what's expected
	python -m pytest test/ | tee $@

## compile lcmodel from source (via github) and copy here
lcmodel/lcmodel:
	git clone https://github.com/schorschinho/LCModel lcmodel-fortran
	make -C lcmodel-fortran
	cp lcmodel-fortran/binaries/linux/lcmodel $@

## GUI program
gui_example:
	./grid.py  -s test/data/siarray.1.1  -r test/data/rorig.nii  -i test/data/roi_pos.txt  --rois roi4 roi6 --gm_mask test/data/gm_sum.nii.gz

test/data/matlab/kspace.1.1: test/genrate_mat.m
	which matlab && \
		matlab -nosplash -nodesktop -r "try;addpath('test'); genrate_mat(); catch e; disp(e); end; quit;" || :

## Standalone programs

# example mkspectrum call
# also see test/data/pos.txt
# x,y 116,124 => 216+1- => 101,93 => sid3: 93,101
# L Caud in /Volumes/Hera/Projects/7TBrainMech/subjs//11743_20190802/slice_PFC/MRSI_roi/13MP20200207/LT
# diff -y <(datamash -t, transpose < out/spectrum.93.101.dir/spreadsheet.csv ) <(datamash -t, transpose < /Volumes/Hera/Projects/7TBrainMech/subjs//11743_20190802/slice_PFC/MRSI_roi/LCModel/v2idxfix/spectrum.93.101.dir/spreadsheet.csv )

out/spectrum.93.101.py: test/data/siarray.1.1 | out
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

out:
	@mkdir -p out

# not needed by track for provenance
test/data/11743_20190802_kspace.1.1: /Volumes/Hera/Projects/7TBrainMech/subjs/11743_20190802/slice_PFC/MRSI_roi/raw//kspace.1.1
	test -r $^ && cp $< $@ || :
