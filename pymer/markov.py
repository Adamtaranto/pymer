# Copyright 2016 Kevin Murray <kdmfoss@gmail.com>
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import print_function, division, absolute_import

import numpy as np
import scipy as sp
from scipy import sparse
from scipy.sparse import linalg

from .base import BaseCounter
from ._hash import (
    hash_to_kmer,
)


class TransitionKmerCounter(BaseCounter):

    '''Counts markovian state transitions in DNA sequences.

    This class counts transtions between (k-1)-mers (or stems) and their
    following bases. This represents the k-1'th order markov process that (may
    have) generated the underlying DNA sequences.

    A normalised, condensed transtion matrix of shape (4^(k-1), 4) or sparse
    complete transtion matrix (shape (4^(k-1), 4^(k-1)) can be returned. In
    addition, the steady-state vector is calculated from the complete transition
    matrix via eigendecomposition.

    Parameters
    ----------
    k : int
        K-mer length
    alphabet : str
        Alphabet over which values are defined, defaults to "ACGT"
    '''

    def __init__(self, k, alphabet="ACGT", array=None, canonical=False):
        BaseCounter.__init__(self, k, alphabet, canonical)
        if array is not None:
            self.array = array
        else:
            self.array = np.zeros((len(alphabet) ** (k-1), len(alphabet)),
                                  dtype=int)
        self.n = 4**(k-1)
        self._transitions = None
        self._P = None

    def _clear(self):
        self._transitions = None
        self._P = None

    @classmethod
    def _kmer2trans(cls, kmer):
        return kmer >> 2, kmer & 0x3

    def _incr(self, kmer, by=1):
        self._clear()
        stem, to = self._kmer2trans(kmer)
        self.array[stem, to] += by

    def _decr(self, kmer, by=1):
        self._clear()
        stem, to = self._kmer2trans(kmer)
        self.array[stem, to] -= by

    @property
    def transitions(self):
        """Dense [k-1]x4 transition frequency matrix"""
        if self._transitions is not None:
            return self._transitions
        transitions = self.array.astype(np.float)
        transitions /= transitions.sum(1)[:, np.newaxis]
        self._transitions = transitions
        return transitions

    @property
    def P(self):
        """Sparse [k-1]x[k-1] transition frequency matrix"""
        if self._P is not None:
            return self._P
        sparse_P = sparse.lil_matrix((self.n, self.n))
        alpha_size = len(self.alphabet)
        bitmask = (self.n-1)  # Mask w/ all bits set hi within (k-1)mer range.
        for fr in range(self.n):
            for a in range(alpha_size):
                to = (fr << 2 | a) & bitmask
                sparse_P[fr, to] = self.transitions[fr, a]
        self._P = sparse_P
        return sparse_P

    @property
    def steady_state(self):
        """Steady-state frequencies of each [k-1]-mer"""
        v, w = linalg.eigs(self.P.transpose(), which='LR')
        ssf = np.real(w[:, v.argmax()])
        ssf /= ssf.sum()
        return ssf

    @property
    def stem_frequencies(self):
        """Frequencies of each stem ([k-1]-mer)"""
        stemfreq = self.array.sum(axis=1).astype(np.float)
        stemfreq /= stemfreq.sum()
        return stemfreq
