cdef class Cache(dict):

    def __getitem__(self, key):
        return dict.__getitem__(self, key)

    def __missing__(self, key):
        ret = self[key] = f(*key)
        return ret

cdef long hash_list(list seq):
    cdef long result = 0x345678
    cdef long mult = 1000003
    cdef long h
    cdef long l = 0
    try:
        l = len(seq)
    except TypeError:
        # NOTE: This probably means very short non-len-able sequences
        # will not be spread as well as they should, but I'm not
        # sure what else to do.
        l = 100
    for element in seq:
        try:
            h = hash(element)
        except TypeError:
            h = hash_list(element)
        result ^= h
        result *= mult
        mult += 82520 + l + l
    result += 97531
    return result ^ hash(type(seq))

