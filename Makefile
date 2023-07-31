.PHONY: test gui_example depends

test: out/test_results.txt

.venv:
	virtualenv .venv/
	source .venv/bin/activate && \
	pip install -r requirements.txt

out/test_results.txt: test/data/11743_20190802_kspace.1.1 test/data/matlab/kspace.1.1 $(wildcard test/*py) lcmodel/lcmodel | out
	# this matlab bit works the same as the original MRRC version
	diff test/data/11743_20190802_kspace.1.1 test/data/11743_20190802_kspace.1.1  | tee -a $@
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
out/spectrum.112.88.py: test/data/siarray.1.1 | out
	@mkdir -p out
	./mkspectrum test/data/siarray.1.1 216 --pos 112 88 out/
	mv out/spectrum.112.88 $@

out/spectrum.112.88.py.dir/: out/spectrum.112.88.py |out
	./lcmodel.py $<

## MATLAB (octave) vesion
out/spectrum.112.88.ml: test/data/siarray.1.1 | out
	octave --eval "addpath('$$PWD/matlab'); gen_spectrum('test/data/siarray.1.1', 216, [112 88],'out/')"
	mv out/spectrum.112.88 $@

out/spectrum.112.88.ml.dir/: out/spectrum.112.88.ml | out
	octave --eval "addpath('$$PWD/matlab/lcmodel/'); cd out/; lcmodel_spec('spectrum.112.88.ml')"

out:
	mkdir -p $@

# not needed by track for provenance
test/data/11743_20190802_kspace.1.1: /Volumes/Hera/Projects/7TBrainMech/subjs/11743_20190802/slice_PFC/MRSI_roi/raw//kspace.1.1
	cp $< $@

