
run_py: interpolate/setup_splines.py
	rm -rf interpolate/__pycache__
	rm -rf interpolate/build
	rm -f **/*.so
	rm -rf **/*.dSYM
	cd interpolate; python setup_splines.py build_ext --inplace
	clear
	clear
	python Spline.py
run_cpp: 