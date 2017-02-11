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
    file_version = 2

    def __init__(self, k, alphabet='ACGT'):
        self.k = k
        self.alphabet = alphabet
        self.num_kmers = len(alphabet) ** k

    def consume(self, seq):
        '''Counts all k-mers in sequence.'''
        for kmer in iter_kmers(seq, self.k):
            self._incr(kmer)

    def consume_file(self, filename):
        """Counts all kmers in all sequences in a FASTA/FASTQ file."""
        with screed.open(filename) as sequences:
            for seq in sequences:
                self.consume(seq['sequence'])

    def unconsume(self, seq):
        '''Subtracts all k-mers in sequence.'''
        for kmer in iter_kmers(seq, self.k):
            self._decr(kmer)

    @classmethod
    def _arraypath(cls, kmersize):
        return "counts_{}-mer".format(kmersize)

    @classmethod
    def read(cls, filename, kmersize):
        h5f = h5py.File(filename, 'r')
        if kmersize not in h5f.attrs['klengths']:
            raise ValueError("Kmer length not in file: {}".format(kmersize))
        arraypath = cls._arraypath(kmersize)
        dset = h5f[arraypath]
        attrs = dset.attrs
        if attrs['class'].decode('utf8') != cls.__name__:
            msg = 'Class mismatch: use {}.read() instead'.format(attrs['class'])
            raise ValueError(msg)
        if attrs['fileversion'] != cls.file_version:
            msg = 'File format version mismatch'
            raise ValueError(msg)
        alphabet = attrs['alphabet'].decode('utf8')
        array = dset[...]
        return cls(kmersize, alphabet=alphabet, array=array)

    def write(self, filename):
        h5f = h5py.File(filename, 'a')
        attrs = {
            'alphabet': self.alphabet,
            'class': self.__class__.__name__,
            'fileversion': self.file_version,
            'pymerversion': get_versions()['version'],
        }
        arraypath = self._arraypath(self.k)
        ds = h5f.create_dataset(arraypath, data=self.array, chunks=True,
                                compression='gzip', compression_opts=9)
        for attr, val in attrs.items():
            if isinstance(val, str):
                val = val.encode('utf-8')
            ds.attrs[attr] = val
        klengths = h5f.attrs.get("klengths", []) + [self.k]
        h5f.attrs['klengths'] = klengths
        h5f.close()

