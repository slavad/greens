target: Makefile
	python setup.py build_ext --inplace

.PHONY: clean

rebuild: clean target

clean:
	rm -f greens/*.so
	rm -rf greens/*.dSYM
	rm -f src/*module.c
	rm -rf src/*wrappers.f
	rm -f  build/*.o
	rm -f  build/*.mod
	rm -rf build/src.*
	rm -rf build/temp.*
	rm -rf build/lib.*

