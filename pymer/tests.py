# Copyright 2016 Kevin Murray <kdmfoss@gmail.com>
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import numpy as np
from tempdir import run_in_tempdir

import itertools as itl

from . import (
    ExactKmerCounter,
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

    for kc in [ExactKmerCounter(2), ]:
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

    for kc in [ExactKmerCounter(3), ]:
        do_test(kc)


@run_in_tempdir()
def test_counter_io():
    for CounterType in [ExactKmerCounter, ]:
        mer = 'AA'

        kc = CounterType(len(mer))

        kc.consume(mer)
        assert kc[mer] == 1

        filename = 'counter.bcz'
        kc.write(filename)
        newkc = CounterType.read(filename, kc.k)
        assert newkc[mer] == 1
