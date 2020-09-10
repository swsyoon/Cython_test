# Import needed setup functions.
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
# We still need to include the directory
# containing the NumPy headers.
#from numpy import get_include
# We still need to run one command via command line.
#from os import system
# Compile the .o file we will be accessing.
# This is independent building the Python extension module.
# shared_obj = "gcc ctridiag.c -fPIC -c -o ctridiag.o"
# print( shared_obj)
# system(shared_obj)
# Tell Python how to compile the extension.
ext_modules = [Extension(
    # Module name:
    name = "cython_ctridiag",
    # Cython source file:
    sources = ["ctridiag.pyx"],
    # Other compile arguments
    libraries=["ctridiag"],
    library_dirs=["lib"],
    include_dirs=["lib"]    
    )]

# Build the extension.
setup(name = 'cython_ctridiag',
      ext_modules = cythonize(ext_modules))


# examples_extension = Extension(
#     name="pyexamples",
#     sources=["pyexamples.pyx"],
#     libraries=["mylib"],
#     library_dirs=["lib"],
#     include_dirs=["lib"]
# )
# setup(
#     name="pyexamples",
#     ext_modules=cythonize(examples_extension)
# )
