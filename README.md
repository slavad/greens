# Intro
Green's function implementation (Fortran) and Python wrapper

# Setup Procedure

This package is currently in the dev stage, so this process will change in the future.

This process is tested on macOS

First of all you must have gfortran installed. If you installed it via `macports` do the following:
  1. `port select --list gcc` to see the list of available gcc versions
  2. `sudo port select --set gcc mp-gcc8` to select `mp-gcc8` for example

Next install Python and `numpy` (this package was tested under 3.6.5).
1. Use `pyenv` to install the version you need: https://github.com/pyenv/pyenv
2. cd to the package dir and run `pip install numpy` to install `numpy`

Now you can build the native extension: cd to the dir of the package and run `make` to build, `make clean` to clean or `make rebuild` to rebuild.


To run example:
  1. ``export PYTHONPATH=`pwd` `` -- run it only once after you open the terminal window. It makes Python to look packages in the current directory.
  2. `python examples/test.py` to run the examples


# Module Structure
  1. `greens.native_functions` contains wrappers for native Fortran functions
  2. `greens.functions` contains additional Python wrappers for functions from `greens.native_functions`.

# File Structure
  `.python-version` contains the Python version used by `pyenv`.

  `src` contains Fortran source code and wrapper declaration (`greens.pyf`)

  `greens` contains Python wrapper functions. Please, note they differ from the wrapper functions defined in `greens.pyf`: sometimes this declaration is
  not enough. E.g. `greens.native_functions.greenlip` takes array `k_i` and it's length `num` as arguments, it's a bit odd.
  So I created a new wrapper `greens.functions.greenlip` which determines the value of `num` and calls `greens.native_functions.greenlip`. Some functions
  are imported from `native_functions` to `functions` without changes. It's recommended to use `greens.functions` for your needs.
