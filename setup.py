# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from os import path

here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

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
    ]
)

from numpy.distutils.core import Extension

ext1 = Extension(name = 'greens',
    sources = [
        'src/GreenFunLip.f90'
    ]
)

from numpy.distutils.core import setup
setup(name = 'py-greens-functions.flib',
    ext_modules = [ext1]
)