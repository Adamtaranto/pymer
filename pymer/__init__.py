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

'''
.. currentmodule:: pymer

This package provides several classes and utilities for counting k-mers in DNA
sequences.

Examples
--------

.. note:: The API demonstrated below applies to all Counters, though Counter
          intialisation varies.

>>> kc = ExactKmerCounter(4)

DNA sequences are counted using the ``consume`` method:

>>> kc.consume('ACGTACGTACGTAC')
>>> kc['ACGT']
3

Sequences can be subtracted using the ``unconsume`` method:

>>> kc.unconsume('ACGTA')
>>> kc['ACGT']
2
>>> kc['CGTA']
2
>>> kc['GTAC']
3

Counters can be added and subtracted:

>>> kc += kc
>>> kc['GTAC']
6
>>> kc -= kc
>>> kc['GTAC']
0

Counters may be read and written to a file, using ``bcolz``.

>>> from tempfile import mkdtemp
>>> from shutil import rmtree
>>> tmpdir = mkdtemp()
>>> filename = tmpdir + '/kc.bcz'

(Above we simply create a temporary directory to hold the saved counts.)

>>> kc.write(filename)
>>> new_kc = ExactKmerCounter.read(filename)
>>> (kc.array == new_kc.array).all()
True
>>> rmtree(tmpdir)


Data Structures
---------------

Summary
^^^^^^^

.. autosummary::

    ExactKmerCounter
    CountMinKmerCounter

Exact K-mer Counting
^^^^^^^^^^^^^^^^^^^^

.. autoclass:: ExactKmerCounter

Probablisistic K-mer Counting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: CountMinKmerCounter

'''

from __future__ import absolute_import, division, print_function
import struct

import bcolz
from bcolz import carray
import numpy as np

from ._hash import (
    iter_kmers,
    hash_to_kmer,
)

from ._cms import (
    cms_getitem,
    cms_setitem,
    cms_incritem,
    cms_decritem,
)

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

__all__ = [
    'ExactKmerCounter',
    'CountMinKmerCounter',
    'iter_kmers',
    'hash_to_kmer',
]


class BaseCounter(object):
    file_version = 1

    def _incr(self, kmer, by=1):
        self[kmer] = min(self.typemax, self[kmer] + by)

    def _decr(self, kmer, by=1):
        self[kmer] = max(self.typemin, self[kmer] - by)

    def consume(self, seq):
        '''Counts all k-mers in sequence.'''
        for kmer in iter_kmers(seq, self.k):
            self._incr(kmer)

    def unconsume(self, seq):
        '''Subtracts all k-mers in sequence.'''
        for kmer in iter_kmers(seq, self.k):
            self._decr(kmer)

    @classmethod
    def read(cls, filename, force_numpy=True):
        array = bcolz.open(filename)
        attrs = array.attrs
        if attrs['class'] != cls.__name__:
            msg = 'Class mismatch: use {}.read() instead'.format(attrs['class'])
            raise ValueError(msg)
        if attrs['fileversion'] != cls.file_version:
            msg = 'File format version mismatch'
            raise ValueError(msg)
        k = attrs['k']
        alphabet = attrs['alphabet']
        if force_numpy and not isinstance(array, np.ndarray):
            array = np.array(array)
        return cls(k, alphabet=alphabet, array=array)

    def write(self, filename):
        attrs = {
            'k': self.k,
            'alphabet': list(self.alphabet),
            'class': self.__class__.__name__,
            'fileversion': self.file_version,
            'pymerversion': __version__,
        }
        array = carray(self.array, rootdir=filename, mode='w')
        for attr, val in attrs.items():
            array.attrs[attr] = val
        array.flush()


class ExactKmerCounter(BaseCounter):

    '''Count k-mers in DNA sequences exactly using an array.

    .. note:: This class is not suitable for k-mers of more than 12 bases.

    Parameters
    ----------
    k : int
        K-mer length
    alphabet : list-like (str, bytes, list, set, tuple) of letters
        Alphabet over which values are defined
    '''

    def __init__(self, k, alphabet='ACGT', array=None):
        self.k = k
        self.alphabet = alphabet
        self.num_kmers = len(alphabet) ** k
        if array is not None:
            self.array = array
        else:
            self.array = np.zeros((len(alphabet) ** k, 1), dtype=int)
        self.typemax = np.iinfo(self.array.dtype).max
        self.typemin = np.iinfo(self.array.dtype).min

    def __add__(self, other):
        if self.k != other.k or self.alphabet != other.alphabet:
            msg = "Cannot add KmerCounters unless k and alphabet are equal."
            raise ValueError(msg)
        x = self.__class__(self.k, self.alphabet)
        x.array = self.array.copy()
        x.array += other.array
        return x

    def __sub__(self, other):
        if self.k != other.k or self.alphabet != other.alphabet:
            msg = "Cannot add KmerCounters unless k and alphabet are equal."
            raise ValueError(msg)
        x = self.__class__(self.k, self.alphabet)
        x.array = self.array.copy()
        x.array -= other.array
        x.array = x.array.clip(min=0)
        return x

    def __len__(self):
        return self.array.sum()

    def __getitem__(self, item):
        if isinstance(item, (str, bytes)):
            if len(item) != self.k:
                msg = "KmerCounter must be queried with k-length kmers"
                return ValueError(msg)
            item = next(iter_kmers(item, self.k))
        return self.array[item, 0]

    def __setitem__(self, item, val):
        if isinstance(item, (str, bytes)):
            if len(item) != self.k:
                msg = "KmerCounter must be queried with k-length kmers"
                return ValueError(msg)
            item = kmer_hash(item)
        self.array[item, 0] = val

    def to_dict(self, sparse=True):
        d = {}
        for kmer in range(self.num_kmers):
            count = self.array[kmer]
            if sparse and count == 0:
                continue
            kmer = hash_to_kmer(kmer, self.k)
            d[kmer] = count
        return d

    def print_table(self, sparse=False, file=None, sep='\t'):
        for kmer, count in sorted(self.to_dict(sparse=sparse).items()):
            print(kmer, count, sep=sep, file=file)


class CountMinKmerCounter(BaseCounter):

    '''
    Count k-mers in DNA sequences using a Count-min Sketch

    Parameters
    ----------
    k : int
        K-mer length
    sketchshape: tuple-like
        Number of tables and table size of the Count-min Sketch. For example,
        sketchshape=(4, 100) makes a Count-min Sketch with 4 tables of 100
        bins.
    alphabet : list-like of letters, optional
        Alphabet over which values are defined. Default ``'ACGT'``.
    dtype: numpy data type
        Count-min Sketch bin data type. Default ``np.uint16``
    '''

    def __init__(self, k, sketchshape=(4, 100000), alphabet='ACGT',
                 dtype=np.uint16, array=None):
        self.k = k
        self.alphabet = alphabet
        self.sketchshape = sketchshape
        if array is not None:
            self.array = array
        else:
            num_tables, table_size = sketchshape
            # We store the sketch transposed to play nicer with bcolz.carray
            self.array = np.zeros((table_size, num_tables), dtype=dtype)
        self.typemax = np.iinfo(self.array.dtype).max
        self.typemin = np.iinfo(self.array.dtype).min

    def __add__(self, other):
        if self.array.shape != other.array.shape or self.k != other.k:
            msg = "Cannot add counters unless k and sketch shape are equal."
            raise ValueError(msg)
        x = self.__class__(self.k, self.sketchshape, self.alphabet,
                           self.array.dtype)
        x.array = self.array.copy()
        dtypemax = np.iinfo(x.array.dtype).max
        overflowidx = (dtypemax - x.array) < other.array
        x.array += other.array
        x.array[overflowidx] = dtypemax
        return x

    def __sub__(self, other):
        if self.array.shape != other.array.shape or self.k != other.k:
            msg = "Cannot add counters unless k and sketch shape are equal."
            raise ValueError(msg)
        x = self.__class__(self.k, self.sketchshape, self.alphabet,
                           self.array.dtype)
        x.array = self.array.copy()
        gtidx = x.array < other.array
        x.array -= other.array
        x.array[gtidx] = 0
        return x

    def __len__(self):
        return self.array.sum(axis=1)[0]

    def __getitem__(self, item):
        if isinstance(item, (str, bytes)):
            if len(item) != self.k:
                msg = "KmerCounter must be queried with k-length kmers"
                return ValueError(msg)
            item = next(iter_kmers(item, self.k))
        return cms_getitem(self.array, item, *self.sketchshape)

    def __setitem__(self, item, value):
        if isinstance(item, (str, bytes)):
            if len(item) != self.k:
                msg = "KmerCounter must be queried with k-length kmers"
                return ValueError(msg)
            item = kmer_hash(item)
        cms_setitem(self.array, item, value, *self.sketchshape)

    def _incr(self, item, by=1):
        #if isinstance(item, (str, bytes)):
        #    if len(item) != self.k:
        #        msg = "KmerCounter must be queried with k-length kmers"
        #        return ValueError(msg)
        #    item = kmer_hash(item)
        cms_incritem(self.array, item, by, *self.sketchshape)

    def _decr(self, item, by=1):
        #if isinstance(item, (str, bytes)):
        #    if len(item) != self.k:
        #        msg = "KmerCounter must be queried with k-length kmers"
        #        return ValueError(msg)
        #    item = kmer_hash(item)
        cms_decritem(self.array, item, by, *self.sketchshape)
