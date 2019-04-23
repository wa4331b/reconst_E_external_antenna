from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os

home_dir = os.environ["HOME"]

ext = Extension(name="poynting",
                sources=["poynting.pyx", "c++/poynting.cpp",
                         "c++/GK_triangle.cpp", "c++/GL.cpp",
                         "c++/triangle_int_FS.cpp"],
                depends=["c++/poynting.hpp", "c++/GK_triangle.hpp",
                         "c++/GL.hpp", "c++/triangle_int.hpp",
                         "c++/dictionary.hpp"],
#                include_dirs=[home_dir+"/includes"],
                include_dirs=["includes"],
                language="c++")
setup(ext_modules=cythonize(ext))
