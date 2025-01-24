"""
Package setup for Cython extensions only.

This script sets up the Cython extensions for the pymer package.

"""

from Cython.Build import build_ext
from setuptools import setup, Extension
import numpy as np

ext_modules = [
    Extension(
        "pymer._hash",
        ["src/pymer/_hash.pyx"],
        include_dirs=[np.get_include()],
        extra_compile_args=["-O3", "-march=native"],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
    ),
]

cmdclass = {"build_ext": build_ext}

setup(
    ext_modules=ext_modules,
    cmdclass=cmdclass,
)