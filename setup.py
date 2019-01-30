# Always prefer setuptools over distutils
from __future__ import division, absolute_import, print_function
from setuptools import setup, find_packages
from os import path
from numpy.distutils.core import Extension
from numpy.distutils.core import setup

here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

ext1 = Extension(name = 'greens',
    sources = [
        'src/GreenFunLip.f90',
        'src/greens.pyf',
        'src/mf_BesselMod1.f90',
        'src/mf_c.f90',
        'src/lib_specfun.f90',
        'src/mf_DIN_djmuz.f'
    ]
)
if __name__ == "__main__":
    setup(
        name='py-greens-functions',
        version='0.0.1',
        description="Python Green's functions implementation",
        long_description=long_description,
        long_description_content_type='text/markdown',
        url='https://github.com/slavad/py-series-clean',
        author='Alexander Mushtukov, Viacheslav Dushin',
        keywords="Green's functions",
        packages=find_packages(),
        classifiers=[
            #TODO: test with other versions
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Fortran :: 77',
            'Programming Language :: Fortran :: 90',
            'Development Status :: 3 - Alpha',
            'License :: OSI Approved :: MIT License'
        ],
        install_requires=[
           'numpy>=1.15.0<2.0.0'
        ],
        ext_modules = [ext1]
    )