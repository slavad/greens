wrapper: Makefile
	python setup.py build_ext --inplace

.PHONY: clean

clean:
	cd build && make clean
	rm -f *.so
	rm -rf *.dSYM

