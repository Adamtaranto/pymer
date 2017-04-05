# Copyright 2016 Kevin Murray <kdmfoss@gmail.com>
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import print_function, division, absolute_import

import numpy as np

from .base import BaseCounter
from ._hash import (
    iter_kmers,
    hash_to_kmer,
)


class ExactKmerCounter(BaseCounter):

    '''Count k-mers in DNA sequences exactly using an array.

    Parameters
    ----------
    k : int
        K-mer length
    alphabet : list-like (str, bytes, list, set, tuple) of letters
        Alphabet over which values are defined, defaults to "ACGT"
    '''

    def __init__(self, k, alphabet='ACGT', array=None, canonical=False):
        BaseCounter.__init__(self, k, alphabet, canonical)
        if array is not None:
            self.array = array
        else:
            self.array = np.zeros(len(alphabet) ** k, dtype=int)
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
            item = next(iter_kmers(item, self.k, self.canonical))
        return self.array[item]

    def __setitem__(self, item, val):
        if isinstance(item, (str, bytes)):
            if len(item) != self.k:
                msg = "KmerCounter must be queried with k-length kmers"
                return ValueError(msg)
            item = next(iter_kmers(item, self.k))
        self.array[item] = val

    def _incr(self, kmer, by=1):
        self[kmer] = min(self.typemax, self[kmer] + by)

    def _decr(self, kmer, by=1):
        self[kmer] = max(self.typemin, self[kmer] - by)

    def to_dict(self, sparse=True):
        d = {}
        for kmer in range(self.num_kmers):
            count = self[kmer]
            if sparse and count == 0:
                continue
            kmer = hash_to_kmer(kmer, self.k)
            d[kmer] = count
        return d

    def print_table(self, sparse=False, file=None, sep='\t'):
        '''
        Prints a table of k-mer counts

        Parameters
        ----------
        sparse: bool
            Print only observed kmers
        file: IO object
            Print to `file`
        sep: str
            Separate k-mers and their counts with `sep`
        '''

        for kmer, count in sorted(self.to_dict(sparse=sparse).items()):
            print(kmer, count, sep=sep, file=file)

    @property
    def counts(self):
        '''k-mer Counts'''
        return self.array

    @property
    def frequencies(self):
        '''k-mer Frequencies'''
        return self.array / self.array.sum()
