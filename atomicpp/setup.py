# Cython compile instructions

from distutils.core import setup, Extension
from Cython.Build import cythonize

# Use python setup.py build_ext --inplace
# to compile

ext_module = Extension(
    "atomicpy",
    ["*.pyx"],
    language="c++",
    extra_compile_args=["-std=c++11"],
    extra_link_args=["-std=c++11"]
)

setup(
  name = "atomicpy",
  ext_modules = cythonize(ext_module),
)