# Copyright 2016 Kevin Murray <kdmfoss@gmail.com>
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

'''
.. currentmodule:: pymer

This package provides several classes and utilities for counting k-mers in DNA
sequences.

Examples
--------

.. note:: The API demonstrated below applies to all Counters, though Counter
          intialisation varies.

>>> ksize = 4
>>> kc = ExactKmerCounter(ksize)

DNA sequences are counted using the ``consume`` method:

>>> kc.consume('ACGTACGTACGTAC')
>>> kc['ACGT']
3

Sequences can be subtracted using the ``unconsume`` method:

>>> kc.unconsume('ACGTA')
>>> kc['ACGT']
2
>>> kc['CGTA']
2
>>> kc['GTAC']
3

Counters can be added and subtracted:

>>> kc += kc
>>> kc['GTAC']
6
>>> kc -= kc
>>> kc['GTAC']
0

Counters may be read and written to a file, using ``HDF5``.

>>> from tempfile import mkdtemp
>>> from shutil import rmtree
>>> tmpdir = mkdtemp()
>>> filename = tmpdir + '/kc.h5'

(Above we simply create a temporary directory to hold the saved counts.)

>>> kc.write(filename)
>>> new_kc = ExactKmerCounter.read(filename, ksize)
>>> (kc.array == new_kc.array).all()
True
>>> rmtree(tmpdir)


Kmers whose strand is unknown (e.g. those from NGS reads) can be "canonicalised"
to the lexographically smaller of the kmer and its reverse complement by
supplying the `canonical=True` argument to `Counter` constructors, or
`iter_kmers`.

>>> canon = ExactKmerCounter(2, canonical=True)
>>> canon.consume('AAAA')
>>> canon['AA']
3
>>> canon['TT']
3
>>> noncanon = ExactKmerCounter(2, canonical=False)
>>> noncanon.consume('AAAA')
>>> noncanon['AA']
3
>>> noncanon['TT']
0


Data Structures
---------------

Summary
^^^^^^^

.. autosummary::

    ExactKmerCounter
    TransitionKmerCounter


Exact K-mer Counting
^^^^^^^^^^^^^^^^^^^^

.. autoclass:: ExactKmerCounter

Markovian K-mer Counting
^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: TransitionKmerCounter

'''

from __future__ import absolute_import, division, print_function

from .base import BaseCounter
from ._hash import (
    iter_kmers,
    hash_to_kmer,
)

from .count import (
    ExactKmerCounter,
)

from .markov import (
    TransitionKmerCounter,
)

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

__all__ = [
    'ExactKmerCounter',
    'TransitionKmerCounter',
    'iter_kmers',
    'hash_to_kmer',
]
