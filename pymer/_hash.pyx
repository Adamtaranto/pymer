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

cdef rjmix(u64 value, u32 seed):
    # An adaptation of Robert Jenkins' 96 bit Mix function
    cdef u64 two32 = 0xFFFFFFFFFF

    # use the upper & lower 32 bits of value as a & b
    cdef u32 a = (value & two32 << 32) >> 32
    cdef u32 b = value & two32
    cdef u32 c = seed

    # below is the original mix function
    a=a-b;  a=a-c;  a=a^(c >> 13);
    b=b-c;  b=b-a;  b=b^(a << 8);
    c=c-a;  c=c-b;  c=c^(b >> 13);
    a=a-b;  a=a-c;  a=a^(c >> 12);
    b=b-c;  b=b-a;  b=b^(a << 16);
    c=c-a;  c=c-b;  c=c^(b >> 5);
    a=a-b;  a=a-c;  a=a^(c >> 3);
    b=b-c;  b=b-a;  b=b^(a << 10);
    c=c-a;  c=c-b;  c=c^(b >> 15);
    return c

def inthash(int value, int seed=0):
    return rjmix(value, seed)
