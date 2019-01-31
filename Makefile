target: Makefile
	python setup.py build_ext --inplace

.PHONY: clean

rebuild: clean target

clean:
	rm -rf greens/native_functions*
	rm -rf src/native_functions*
	rm -rf build/*

