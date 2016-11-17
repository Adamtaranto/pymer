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

Counters may be read and written to a file, using ``bcolz``.

>>> from tempfile import mkdtemp
>>> from shutil import rmtree
>>> tmpdir = mkdtemp()
>>> filename = tmpdir + '/kc.bcz'

(Above we simply create a temporary directory to hold the saved counts.)

>>> kc.write(filename)
>>> new_kc = ExactKmerCounter.read(filename, ksize)
>>> (kc.array == new_kc.array).all()
True
>>> rmtree(tmpdir)


Data Structures
---------------

Summary
^^^^^^^

.. autosummary::

    ExactKmerCounter

Exact K-mer Counting
^^^^^^^^^^^^^^^^^^^^

.. autoclass:: ExactKmerCounter

'''

from __future__ import absolute_import, division, print_function

from ._hash import (
    iter_kmers,
    hash_to_kmer,
)

from .count import (
    ExactKmerCounter,
)

from .markov import (
    KmerTransitionCounter,
)

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

__all__ = [
    'ExactKmerCounter',
    'iter_kmers',
    'hash_to_kmer',
]
