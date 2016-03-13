'''
Nucleotide word (k-mer) counters (:mod:`pymer`)
===============================================

.. currentmodule:: pymer

This package provides several classes and utilities for counting k-mers in DNA
sequences.

Data Structures
---------------

.. autosummary::
    :toctree: generated/

    ExactKmerCounter
    CountMinKmerCounter

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

Counters may be read and written to a file. A YAML object describing the
counter is writen, with arrays compressed with bloscpack.

>>> dumped = kc.write()
>>> new_kc = ExactKmerCounter.read(string=dumped)
>>> (kc.array == new_kc.array).all()
True
'''

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

from __future__ import absolute_import, division, print_function
import yaml
import struct

import bloscpack as bp
import numpy as np
from xxhash import xxh64

from ._hash import (
    iter_kmers,
    hash_to_kmer,
)


__all__ = [
    'ExactKmerCounter',
    'CountMinKmerCounter',
    'iter_kmers',
    'hash_to_kmer',
]


class BaseCounter(object):
    file_version = 1

    def consume(self, seq):
        '''Counts all k-mers in sequence.'''
        for kmer in iter_kmers(seq, self.k):
            self[kmer] += 1

    def unconsume(self, seq):
        '''Subtracts all k-mers in sequence.'''
        for kmer in iter_kmers(seq, self.k):
            self[kmer] = max(self[kmer] - 1, 0)

    @classmethod
    def read(cls, filename=None, string=None):

        if filename is not None:
            with open(filename) as fh:
                obj = yaml.load(fh)
        else:
            obj = yaml.load(string)
        if obj['class'] != cls.__name__:
            msg = 'Class mismatch: use {}.read() instead'.format(obj['class'])
            raise ValueError(msg)
        if obj['fileversion'] < cls.file_version:
            msg = 'File format version mismatch'
            raise ValueError(msg)
        k = obj['k']
        alphabet = obj['alphabet']
        array = bp.unpack_ndarray_str(obj['array'])
        return cls(k, alphabet=alphabet, array=array)


    def write(self, filename=None):
        array = bp.pack_ndarray_str(self.array)
        obj = {'k': self.k,
               'alphabet': list(self.alphabet),
               'array': array,
               'class': self.__class__.__name__,
               'fileversion': self.file_version,
               }

        if filename is not None:
            with open(filename, 'w') as fh:
                yaml.dump(obj, fh)
        else:
            return yaml.dump(obj)





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
    writables = ['k', 'alphabet', 'array']

    def __init__(self, k, alphabet='ACGT', array=None):
        self.k = k
        self.alphabet = alphabet
        self.num_kmers = len(alphabet) ** k
        if array is not None:
            self.array = array
        else:
            self.array = np.zeros(len(alphabet) ** k, dtype=int)

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
        return self.array[item]

    def __setitem__(self, item, val):
        if isinstance(item, (str, bytes)):
            if len(item) != self.k:
                msg = "KmerCounter must be queried with k-length kmers"
                return ValueError(msg)
            item = kmer_hash(item)
        self.array[item] = val

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
        self.num_tables, self.table_size = sketchshape
        if array is not None:
            self.array = array
        else:
            self.array = np.zeros(sketchshape, dtype=dtype)

    def __add__(self, other):
        if self.array.shape != other.array.shape or self.k != other.k:
            msg = "Cannot add counters unless k and sketch shape are equal."
            raise ValueError(msg)
        x = self.__class__(self.k, self.array.shape, self.alphabet,
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
        x = self.__class__(self.k, self.array.shape, self.alphabet,
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
        mx = 0
        for tab in range(self.num_tables):
            idx = xxh64(struct.pack('Q', item), seed=tab).intdigest()
            idx %= self.table_size
            mx = max(mx, self.array[tab, idx])
        return mx

    def __setitem__(self, item, val):
        if isinstance(item, (str, bytes)):
            if len(item) != self.k:
                msg = "KmerCounter must be queried with k-length kmers"
                return ValueError(msg)
            item = kmer_hash(item)
        for tab in range(self.num_tables):
            idx = xxh64(struct.pack('Q', item), seed=tab).intdigest()
            idx %= self.table_size
            self.array[tab, idx] = val
