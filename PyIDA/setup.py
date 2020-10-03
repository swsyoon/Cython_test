from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

ext_modules=[   
    Extension("sundials",
              sources = ['IDA.pyx'],
              include_dirs = ["include"],
              library_dirs=["lib"],
              libraries = ['sundials_ida','sundials_nvecserial'])]
              

# Build the extension.
setup(name = 'cython_IDA', ext_modules = cythonize(ext_modules))



#!python setup.py build_ext --inplace