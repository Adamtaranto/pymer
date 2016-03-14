import struct
import numpy as np
cimport numpy as np


ctypedef unsigned short u16
ctypedef u16 DTYPE
ctypedef unsigned long long int u64

cdef extern from "xxhash.h":
    u64 XXH64 (const void* input, size_t length, u64 seed)


cdef inline u64 u64hash(u64 val, u64 seed):
    # This is the xxhash of a u64
    cdef unsigned char *b = <unsigned char *>&val
    return XXH64(b, 8, seed)


def cms_getitem(np.ndarray[DTYPE, ndim=2] array not None, u64 item, u64 ntab, u64 tabsize):
    cdef u64 mx = 0
    cdef u64 idx, tab
    cdef u16 array_val
    for tab in range(ntab):
        idx = u64hash(item, tab) % tabsize
        array_val = array[tab, idx]
        if mx < array_val:
            mx = array_val
    return mx


def cms_setitem(np.ndarray[DTYPE, ndim=2] array not None, u64 item, u64 value, u64 ntab,
                u64 tabsize):
    cdef u64 idx, tab
    for tab in range(ntab):
        idx = u64hash(item, tab) % tabsize
        array[tab, idx] = value

