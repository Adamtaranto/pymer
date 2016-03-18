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
cimport numpy as np
cimport cython


ctypedef unsigned short u16
ctypedef u16 DTYPE
ctypedef unsigned long long int u64

cdef extern from "xxhash.h":
    u64 XXH64 (const void* input, size_t length, u64 seed)


cdef inline u64 u64hash(u64 val, u64 seed):
    # This is the xxhash of a u64
    cdef unsigned char *b = <unsigned char *>&val
    return XXH64(b, 8, seed)


def cms_incritem(np.ndarray[DTYPE, ndim=2] array not None, u64 item, DTYPE by,
                 u64 ntab, u64 tabsize):
    cdef u64 idx, tab
    cdef DTYPE array_val, orig_val
    for tab in range(ntab):
        idx = u64hash(item, tab) % tabsize
        orig_val = array[idx, tab]
        array[idx, tab] += by
        if array[idx, tab] < orig_val:
            array[idx, tab] = orig_val

def cms_decritem(np.ndarray[DTYPE, ndim=2] array not None, u64 item, DTYPE by,
                 u64 ntab, u64 tabsize):
    cdef u64 idx, tab
    cdef DTYPE array_val, orig_val
    for tab in range(ntab):
        idx = u64hash(item, tab) % tabsize
        orig_val = array[idx, tab]
        array[idx, tab] -= by
        if array[idx, tab] > orig_val:
            array[idx, tab] = orig_val

def cms_getitem(np.ndarray[DTYPE, ndim=2] array not None, u64 item, u64 ntab,
                u64 tabsize):
    cdef u64 mx = 0
    cdef u64 idx, tab
    cdef DTYPE array_val
    for tab in range(ntab):
        idx = u64hash(item, tab) % tabsize
        array_val = array[idx, tab]
        if mx < array_val:
            mx = array_val
    return mx


def cms_setitem(np.ndarray[DTYPE, ndim=2] array not None, u64 item,
                DTYPE value, u64 ntab, u64 tabsize):
    cdef u64 idx, tab
    for tab in range(ntab):
        idx = u64hash(item, tab) % tabsize
        array[idx, tab] = value

