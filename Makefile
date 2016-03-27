
.PHONY: all test
all:
	python setup.py build_ext --inplace

test:
	./setup.py test

htmlcov:
	nosetests --with-cov --cov-report html --cover-package pymer
