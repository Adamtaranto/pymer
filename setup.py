#!/usr/bin/env python3
#
# Copyright 2016 Kevin Murray <spam@kdmurray.id.au>
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from setuptools import setup
import versioneer
from setuptools.extension import Extension
from setuptools.command.test import test as TestCommand
from setuptools.command.build_ext import build_ext as OrigBuildExt

EXT = 'c'
try:
    from Cython.Build import cythonize
    EXT = 'pyx'
except ImportError:
    def cythonize(x): return x

description = """
pymer: Pythonic fast k-mer counting routines
"""


class NoseCommand(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        # Run nose ensuring that argv simulates running nosetests directly
        import nose
        nose.run_exit(argv=['nosetests'])


class BuildExt(OrigBuildExt):

    def finalize_options(self):
        OrigBuildExt.finalize_options(self)

        import numpy
        self.include_dirs.append(numpy.get_include())


cmdclasses = versioneer.get_cmdclass()
cmdclasses['test'] = NoseCommand
cmdclasses['build_ext'] = BuildExt

setup(
    name='pymer',
    packages=['pymer', ],
    version=versioneer.get_version(),
    cmdclass=cmdclasses,
    install_requires=[
        'bcolz',
    ],
    setup_requires=[
        'cython>=0.23',
        'numpy>=1.8',
        'nose',
    ],
    tests_require=[
        'nose-cov',
        'blessings',
        'docopt',
        'tempdir',
    ],
    ext_modules=cythonize([
        Extension(
            'pymer._hash', [
                'pymer/_hash.{}'.format(EXT),
            ],
        ),
        Extension(
            'pymer._cms', [
                'pymer/_cms.{}'.format(EXT),
                'pymer/xxhash.c',
             ],
        ),
    ]),
    description=description,
    author="Kevin Murray",
    author_email="spam@kdmurray.id.au",
    url="https://github.com/kdmurray91/pymer",
    keywords=["kmer"],
    classifiers=[
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
