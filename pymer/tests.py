import numpy as np
import itertools as itl

from . import (
    KmerCounter,
)
from ._pymer import (
    iter_kmers,
    hash_to_kmer,
)


# de Bruijn DNA sequences of k={2,3}, i.e. contain all 2/3-mers once
K2_DBS = 'AACAGATCCGCTGGTTA'
K3_DBS = 'AAACAAGAATACCACGACTAGCAGGAGTATCATGATTCCCGCCTCGGCGTCTGCTTGGGTGTTTAA'

def test_counter_init():
    kc = KmerCounter(5)
    assert kc.k == 5
    assert kc.num_kmers == 4**5
    assert list(kc.alphabet) == list('ACGT')
    assert np.all(kc.array == np.zeros(4**5, dtype=int))
    assert len(kc) == 0

    kc = KmerCounter(5, alphabet='NOTDNA')
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


def test_hash_to_kmer():
    k = 2
    hashes = range(4**k)
    kmers = map(''.join, list(itl.product(list('ACGT'), repeat=k)))
    for hsh, mer in zip(hashes, kmers):
        h2k = hash_to_kmer(hsh, k)
        assert h2k == mer, (hsh, mer, h2k)


def test_counter_operations():
    kc = KmerCounter(2)
    kc.consume(K2_DBS)

    add = kc + kc
    assert np.all(add.array == np.ones(4**2, dtype=int)*2)

    sub = add - kc
    assert np.all(sub.array == kc.array)

    sub -= kc
    sub -= kc
    assert np.all(sub.array == np.zeros(4**2, dtype=int))


def test_counter_consume():
    # de Bruijn DNA sequence of k=3, i.e. contains all 3-mers once
    kc = KmerCounter(3)
    kc.consume(K3_DBS)
    assert np.all(kc.array == np.ones(4**3, dtype=int))

    kc.unconsume('ACT')
    assert kc['ACT'] == 0

