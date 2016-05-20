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

import numpy as np
from tempdir import run_in_tempdir

import itertools as itl

from . import (
    ExactKmerCounter,
    CountMinKmerCounter,
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

    for kc in [ExactKmerCounter(2), CountMinKmerCounter(2, (4, 100000))]:
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

    for kc in [ExactKmerCounter(3), CountMinKmerCounter(3, (4, 100000))]:
        do_test(kc)


def test_cms_counter_overflow():
    K = 2
    kc = CountMinKmerCounter(K, (2, 100))

    for _ in range(2**16):
        kc.consume('AA')

    assert kc['AA'] == 2**16 - 1, kc['AA']


@run_in_tempdir()
def test_counter_io():
    for CounterType in ExactKmerCounter, CountMinKmerCounter:
        mer = 'AA'

        kc = CounterType(len(mer))

        kc.consume(mer)
        assert kc[mer] == 1

        filename = 'counter.bcz'
        kc.write(filename)
        newkc = CounterType.read(filename)
        assert kc[mer] == 1
