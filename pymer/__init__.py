from __future__ import absolute_import, division, print_function
import numpy as np

from ._pymer import (
    iter_kmers,
    hash_to_kmer,
)

__all__ = [
    'KmerCounter',
    'iter_kmers',
    'hash_to_kmer',
]

class KmerCounter(object):
    '''
    Count k-mers in DNA sequences.

    DNA sequences are counted using the .consume method:

    >>> kc = KmerCounter(4)
    >>> kc.consume('ACGTACGTACGT')
    >>> kc['ACGT']
    3

    KmerCounters can be added and subtracted, which 
    >>> kc += kc
    >>> kc['ACGT']
    6
    >>> kc -= kc
    >>> kc['ACGT']
    0

    '''

    def __init__(self, k, alphabet='ACGT'):
        self.k = k
        self.alphabet = alphabet
        self.num_kmers = len(alphabet) ** k
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


    def consume(self, seq):
        '''Counts all k-mers in sequence.

        >>> kc = KmerCounter(2)
        >>> kc.consume('AAAA')
        >>> kc['AA']
        3
        '''
        for kmer in iter_kmers(seq, self.k):
            self.array[kmer] += 1

    def unconsume(self, seq):
        '''Subtracts all k-mers in sequence.
        >>> kc = KmerCounter(2)
        >>> kc.consume('AAAA')
        >>> kc['AA']
        3
        >>> kc.unconsume('AA')
        >>> kc['AA']
        2

        Never gives negative numbers:
        >>> kc['AA']
        2
        >>> kc.unconsume('AAAA')
        >>> kc['AA']
        0
        '''
        for kmer in iter_kmers(seq, self.k):
            self.array[kmer] = max(self.array[kmer] - 1, 0)

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
