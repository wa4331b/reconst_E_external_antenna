from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os

home_dir = os.environ["HOME"]

ext = Extension(name="Z_obs",
                sources=["Z_obs.pyx", "c++/compute_Z_obs.cpp",
                         "c++/triangle_int_FS.cpp", "c++/GK_triangle.cpp", 
                         "c++/GL.cpp"],
                depends=["c++/EMConstants.hpp", "c++/GK_triangle.hpp",
                         "c++/GL.hpp", "c++/triangle_int.hpp",
                         "c++/compute_Z_obs.hpp"],
#                include_dirs=[home_dir + "/includes/"],
                include_dirs=["includes/"],
                language="c++")
setup(ext_modules=cythonize(ext))
