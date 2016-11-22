# Copyright 2016 Kevin Murray <kdmfoss@gmail.com>
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import print_function, division, absolute_import

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

    pass
