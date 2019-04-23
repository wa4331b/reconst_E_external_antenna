from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os

home_dir = os.environ["HOME"]

ext = Extension(name="make_trans",
                sources=["make_trans.pyx", "c++/make_trans.cpp", "c++/interpolation.cpp"],
                depends=["c++/make_trans.hpp", "c++/interpolation.hpp", "setup_trans.py"],
#                include_dirs=[home_dir + "/includes/"],
                include_dirs=["includes/"],
                language="c++")
setup(ext_modules=cythonize(ext))
