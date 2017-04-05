# -*- coding: utf-8 -*-
# Copyright 2016 Kevin Murray <kdmfoss@gmail.com>
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


import numpy as np
from tempdir import run_in_tempdir

import itertools as itl

from . import (
    ExactKmerCounter,
    TransitionKmerCounter
)
from ._hash import (
    iter_kmers,
    hash_to_kmer,
)


# de Bruijn DNA sequences of k={2,3}, i.e. contain all 2/3-mers once
K2_DBS = 'AACAGATCCGCTGGTTA'
K3_DBS = 'AAACAAGAATACCACGACTAGCAGGAGTATCATGATTCCCGCCTCGGCGTCTGCTTGGGTGTTTAA'


def all_kmers(k):
    for kmer in itl.product('ACGT', repeat=k):
        yield ''.join(kmer)


def test_counter_init():
    kc = ExactKmerCounter(5)
    assert kc.k == 5
    assert kc.num_kmers == 4**5
    assert list(kc.alphabet) == list('ACGT')
    assert np.all(kc.array == np.zeros(4**5, dtype=int))
    assert len(kc) == 0

    kc = ExactKmerCounter(5, alphabet='NOTDNA')
    assert kc.k == 5
    assert kc.num_kmers == 6**5
    assert list(kc.alphabet) == list('NOTDNA')
    assert np.all(kc.array == np.zeros(6**5, dtype=int))


def test_iter_kmers():
    k = 2
    counts = np.zeros(4**k, dtype=int)
    for kmer in iter_kmers(K2_DBS, k):
        counts[kmer] += 1
    assert counts.sum() == len(K2_DBS) - k + 1, counts.sum()
    assert (counts == 1).all(), counts


def test_iter_kmers_canonical():
    seq = "AAAAAAA"
    rcseq = "TTTTTTT"
    k = 2

    canon_seq = set(iter_kmers(seq, k, canonical=True))
    canon_rcseq = set(iter_kmers(rcseq, k, canonical=True))
    print(canon_seq, canon_rcseq)
    assert canon_seq == canon_rcseq

    non_canon_seq = set(iter_kmers(seq, k, canonical=False))
    non_canon_rcseq = set(iter_kmers(rcseq, k, canonical=False))
    assert non_canon_seq != non_canon_rcseq


def test_iter_kmers_ns():
    k = 3
    seq = "ACGTNACGTNCG"
    expect = [0b000110, 0b011011, 0b000110, 0b011011, ]
    got = list(iter_kmers(seq, k))
    assert got == expect, (got, expect)


def test_hash_to_kmer():
    k = 2
    hashes = range(4**k)
    kmers = map(''.join, list(itl.product(list('ACGT'), repeat=k)))
    for hsh, mer in zip(hashes, kmers):
        h2k = hash_to_kmer(hsh, k)
        assert h2k == mer, (hsh, mer, h2k)


def test_counter_operations():
    def do_test(kc):
        kc.consume(K2_DBS)

        for mer in all_kmers(2):
            assert kc[mer] == 1

        add = kc + kc
        for mer in all_kmers(2):
            assert add[mer] == 2  # each kmer twice

        sub = add - kc
        for mer in all_kmers(2):
            assert sub[mer] == 1  # back to once

        sub -= kc
        sub -= kc
        for mer in all_kmers(2):
            assert sub[mer] == 0, (sub[mer], kc)  # caps at zero even after -2

    for kc in [ExactKmerCounter(2, canonical=False), ]:
        do_test(kc)


def test_counter_consume():
    def do_test(kc):
        for mer in all_kmers(3):
            assert kc[mer] == 0  # zero at start

        kc.consume(K3_DBS)
        for mer in all_kmers(3):
            assert kc[mer] == 1  # After consuming

        kc.unconsume(K3_DBS)
        for mer in all_kmers(3):
            assert kc[mer] == 0  # back to zero after unconsume

    for kc in [ExactKmerCounter(3, canonical=False), ]:
        do_test(kc)


@run_in_tempdir()
def test_counter_readall():
    filename = "counter.h5"
    for CounterType in [ExactKmerCounter, TransitionKmerCounter]:
        for ksize in range(1, 4):
            kc = CounterType(ksize)

            kc.consume(K3_DBS)
            assert kc.array.sum() == len(K3_DBS) - ksize + 1

            kc.write(filename)

        allctrs = CounterType.readall(filename)
        for ksize in range(1, 4):
            assert ksize in allctrs
            assert allctrs[ksize].array.sum() == len(K3_DBS) - ksize + 1


@run_in_tempdir()
def test_counter_io():
    for CounterType in [ExactKmerCounter, TransitionKmerCounter]:
        mer = 'AA'

        kc = CounterType(len(mer))

        kc.consume(mer)
        assert kc.array.sum() == 1

        filename = 'counter.h5'
        kc.write(filename)
        newkc = CounterType.read(filename, kc.k)
        assert newkc.array.sum() == 1


def test_transition_counter():
    '''TransitionKmerCounter behaves correctly'''
    t = TransitionKmerCounter(3, canonical=False)
    t.consume(K3_DBS)

    # Raw counts
    counts = t.array
    assert (counts == 1).all(), counts

    # P, i.e. (k-1, k-1)-sized array of probabilies, the formal transition
    # matrix
    P = t.transitions
    expect_P = np.zeros_like(P) + 0.25  # should all be eq. emission prob
    assert np.allclose(P, expect_P)
    assert (P.sum(1) == 1).all(), P.sum(1)

    # Steady-state frequency vector, π
    pi = t.steady_state
    assert pi.sum() == 1, pi.sum()
    assert np.allclose(pi.dot(t.P.toarray()),  pi), pi

    # Stem frequencies
    sf = t.stem_frequencies
    assert np.allclose(sf, np.ones_like(sf) / len(sf))


def test_exact_counter():
    '''ExactKmerCounter() behaves correctly'''
    k = 3
    c = ExactKmerCounter(k, canonical=False)
    c.consume(K3_DBS)

    assert c.array.shape == (4**k,)

    assert c.array.sum() == len(K3_DBS) - k + 1
    assert len(c) == len(K3_DBS) - k + 1

    assert all(c.counts == 1)
    assert all(c.frequencies == 1/(4**k))
