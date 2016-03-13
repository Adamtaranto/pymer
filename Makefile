
.PHONY: all test2 test3 test
all:
	python setup.py build_ext --inplace
	python3 setup.py build_ext --inplace

test2: all
	nosetests

test3: all
	nosetests3

test: test2 test3
