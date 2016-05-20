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
cimport numpy as cnp
cimport cython


ctypedef unsigned long long int u64
ctypedef unsigned int u32


@cython.boundscheck(False)
def hash_to_kmer(int h, int k):
    '''Convert a hash value at given k to the string representation.
    '''
    cdef const char *nts = 'ACGT'
    cdef u64 nt
    kmer = []
    for x in range(k):
        nt = (h >> (2*x)) & 0x03
        kmer.append(chr(nts[nt]))
    return ''.join(reversed(kmer))


def iter_kmers(str seq, int k):
    '''Iterator over hashed k-mers in a string DNA sequence.
    '''
    cdef u64 n
    cdef u64 bitmask = 2**(2*k)-1  # Set lowest 2*k bits
    cdef u64 h = 0

    # For each kmer's end nucleotide, bit-shift, add the end and yield
    cdef u64 skip = 0
    for end in range(len(seq)):
        nt = seq[end]
        if skip > 0:
            skip -= 1
        if nt == 'A' or nt == 'a':
            n = 0
        elif nt == 'C' or nt == 'c':
            n = 1
        elif nt == 'G' or nt == 'g':
            n = 2
        elif nt == 'T' or nt == 't':
            n = 3
        else:
            skip = k
            continue
        h = ((h << 2) | n) & bitmask
        if end >= k - 1 and skip == 0:
            # Only yield once an entire kmer has been loaded into h
            yield h
