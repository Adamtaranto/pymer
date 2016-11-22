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
    pass
