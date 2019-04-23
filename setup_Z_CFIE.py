from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os

home_dir = os.environ["HOME"]

ext = Extension(name="Z_CFIE",
                sources=["Z_CFIE.pyx", "c++/triangle_int_FS.cpp",
                         "c++/GK_triangle.cpp", "c++/GL.cpp",
                         "c++/Z_EJ_Z_HJ_FS_triangles_arrays.cpp"],
                depends=["dictionary.hpp", "EMConstants.hpp",
                         "GK_triangle.hpp", "GL.hpp",
                         "triangle_int.hpp", "Z_EJ_Z_HJ.hpp"],
#                include_dirs=[home_dir+"/includes/"],
                include_dirs=["includes/"],
                # extra_compile_args=["-std=c++11"],
                language="c++")
setup(ext_modules=cythonize(ext))
