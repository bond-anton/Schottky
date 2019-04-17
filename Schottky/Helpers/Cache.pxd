cdef class Cache(dict):
    cdef:
        unsigned long __hits
        unsigned long __misses

cdef long hash_list(list seq)
