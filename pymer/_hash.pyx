# cython: profile=False

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
    for end in range(len(seq)):
        n = (ord(seq[end]) & 6) >> 1
        n ^= n>>1
        h = ((h << 2) | n) & bitmask
        if end >= k - 1:
            # Only yield once an entire kmer has been loaded into h
            yield h
