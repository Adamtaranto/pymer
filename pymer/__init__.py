from __future__ import absolute_import, division, print_function
import numpy as np

from ._pymer import (
    iter_kmers,
)

class KmerCounter(object):

    def __init__(self, k, alphabet='ACGT'):
        self.k = k
        self.alphabet = alphabet
        self.num_kmers = len(alphabet) ** k
        self.array = np.zeros(len(alphabet) ** k)

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


    def consume(self, seq):
        for kmer in iter_kmers(seq, self.k):
            self.array[kmer] += 1

    def unconsume(self, seq):
        for kmer in iter_kmers(seq):
            self.array[item] -= 1

    def to_dict(self, sparse=True):
        d = {}
        for kmer in range(self.num_kmers):
            count = self.array[kmer]
            if sparse and count == 0:
                continue
            kmer = hash_to_kmer(kmer, self.k)
            d[kmer] = count

    def print_table(self, sparse=False, file=None, sep='\t'):
        for kmer, count in self.to_dict(sparse=sparse).items():
            print(kmer, count, sep=sep, file=file)
