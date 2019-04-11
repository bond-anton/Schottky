cdef class Cache(dict):
    cdef:
        unsigned long __max_size
        unsigned long __id
        unsigned long __ttl
        dict __records, __reverse_records, __records_time
        unsigned long __hits
        unsigned long __misses
        unsigned long __expired

cdef long hash_list(list seq)
