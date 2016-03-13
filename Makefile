
.PHONY: all test2 test3 test
all:
	./setup.py build_ext --inplace

test2:
	nosetests

test3:
	nosetests3

test: test2 test3
