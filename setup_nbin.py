from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os

home_dir = os.environ["HOME"]

ext = Extension(name="nbin",
                sources=["nbin.pyx", "c++/nbin.cpp"],
                depends=["nbin.hpp"],
#                include_dirs=[home_dir+"/includes"],
                include_dirs=["includes"],
                language="c++")
setup(ext_modules=cythonize(ext))
