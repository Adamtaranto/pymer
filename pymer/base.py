# Copyright 2016 Kevin Murray <kdmfoss@gmail.com>
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import numpy as np
import h5py
import screed

from ._hash import (
    iter_kmers,
    hash_to_kmer,
)

from ._version import get_versions


class BaseCounter(object):
    file_version = 4

    def __init__(self, k, alphabet='ACGT', canonical=False):
        self.k = k
        self.alphabet = alphabet
        self.num_kmers = len(alphabet) ** k
        self.canonical = canonical

    def consume(self, seq):
        '''Counts all k-mers in sequence.'''
        for kmer in iter_kmers(seq, self.k, canonical=self.canonical):
            self._incr(kmer)

    def consume_file(self, filename):
        """Counts all kmers in all sequences in a FASTA/FASTQ file."""
        with screed.open(filename) as sequences:
            for seq in sequences:
                self.consume(seq['sequence'])

    def unconsume(self, seq):
        '''Subtracts all k-mers in sequence.'''
        for kmer in iter_kmers(seq, self.k, canonical=self.canonical):
            self._decr(kmer)

    @classmethod
    def _arraypath(cls, kmersize):
        return "{}_{}".format(cls.__name__, kmersize)

    @classmethod
    def read(cls, filename, kmersize):
        h5f = h5py.File(filename, 'r')
        if h5f.attrs['fileversion'] != cls.file_version:
            msg = 'File format version mismatch'
            raise ValueError(msg)
        arraypath = cls._arraypath(kmersize)
        if arraypath not in h5f:
            raise KeyError("k-mer size and counter type not in file")
        dset = h5f[arraypath]
        alphabet = dset.attrs['alphabet'].decode('utf8')
        array = dset[...]
        return cls(kmersize, alphabet=alphabet, array=array)

    @classmethod
    def readall(cls, filename):
        h5f = h5py.File(filename, 'r')
        instances = {}
        for arraypath in h5f.keys():
            clsname, _, k = arraypath.partition('_')
            kmersize = int(k)
            if clsname == cls.__name__:
                instances[kmersize] = cls.read(filename, kmersize)
        return instances

    def write(self, filename):
        h5f = h5py.File(filename, 'a')
        if 'fileversion' in h5f.attrs and \
                h5f.attrs['fileversion'] != self.file_version:
            msg = 'File format version mismatch'
            raise ValueError(msg)

        h5f.attrs['fileversion'] = self.file_version
        h5f.attrs['pymerversion'] = get_versions()['version'].encode('utf8')

        ds = h5f.create_dataset(self._arraypath(self.k), data=self.array,
                                chunks=True, compression='gzip',
                                compression_opts=9)

        ds.attrs['alphabet'] = self.alphabet.encode('utf8')
        h5f.close()
