# Cython compile instructions

from distutils.core import setup, Extension
from Cython.Build import cythonize

# Use python setup.py build_ext --inplace
# to compile

ext_module = Extension(
    "atomicpy",
    ["atomicpy.pyx"],
    language="c++",
    extra_compile_args=['-std=c++11','-O0'],
	extra_link_args=['-std=c++11','-O0']
    # extra_compile_args=['-std=c++11', '-fno-inline', '-g', '-Wall', '-O0'],
    # extra_link_args=['-std=c++11', '-fno-inline', '-g', '-Wall', '-O0']
)

setup(
  name = "atomicpy",
  ext_modules = cythonize(ext_module),
)