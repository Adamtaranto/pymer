from xxhash import xxh64
import struct

def cms_getitem(array, item):
    mx = 0
    ntab, tabsize = array.shape
    for tab in range(ntab):
        idx = xxh64(struct.pack('Q', item), seed=tab).intdigest()
        idx %= tabsize
        mx = max(mx, array[tab, idx])
    return mx

def cms_setitem(array, item, value):
    ntab, tabsize = array.shape
    for tab in range(ntab):
        idx = xxh64(struct.pack('Q', item), seed=tab).intdigest()
        idx %= tabsize
        array[tab, idx] = value

