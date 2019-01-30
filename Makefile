target: Makefile
	python setup.py build_ext --inplace

.PHONY: clean

clean:
	rm -f *.so
	rm -rf *.dSYM
	rm -f src/*module.c
	rm -rf src/*wrappers.f
	rm -f  build/*.o
	rm -f  build/*.mod
	rm -rf build/src.*
	rm -rf build/temp.*
	rm -rf build/lib.*

