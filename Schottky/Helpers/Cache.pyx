from libc.time cimport time

cdef class Cache(dict):

    def __init__(self, max_size=0, ttl=0, *args, **kwargs):
        super(Cache, self).__init__(*args, **kwargs)
        self.__max_size = <unsigned long> abs(max_size)
        self.__id = 0
        self.__ttl = <unsigned long> ttl
        self.__records = {}
        self.__reverse_records = {}
        self.__records_time = {}
        self.__hits = 0
        self.__misses = 0
        self.__expired = 0

    def __getitem__(self, key):
        self.__hits += 1
        return super(Cache, self).__getitem__(key)

    def __missing__(self, key):
        self.__misses += 1
        raise KeyError(key)

    def __setitem__(self, key, value):
        self.__records[key] = self.__id
        self.__reverse_records[self.__id] = key
        self.__records_time[key] = <unsigned long> time(NULL)
        super(Cache, self).__setitem__(key, value)
        if len(self) == self.__max_size + 1:
            self.__delitem__()

    def update(self, *args, **kwargs):
        if args:
            if len(args) > 1:
                raise TypeError("update expected at most 1 arguments, "
                                "got %d" % len(args))
            other = dict(args[0])
            for key in other:
                self[key] = other[key]
        for key in kwargs:
            self[key] = kwargs[key]

    def setdefault(self, key, value=None):
        if key not in self:
            self[key] = value
        return self[key]

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

