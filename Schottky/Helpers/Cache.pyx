cdef class Cache(dict):

    def __init__(self, label=None, *args, **kwargs):
        super(Cache, self).__init__(*args, **kwargs)
        self.__hits = 0
        self.__misses = 0
        if label is None:
            self.__label = 'Cache'
        else:
            self.__label = str(label)

    @property
    def label(self):
        return self.__label

    @label.setter
    def label(self, str label):
        self.__label = label

    def __getitem__(self, key):
        self.__hits += 1
        return super(Cache, self).__getitem__(key)

    def __missing__(self, key):
        self.__misses += 1
        raise KeyError(key)

    def __setitem__(self, key, value):
        super(Cache, self).__setitem__(key, value)

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

    def info(self):
        print('%s:' % self.__label)
        print('  Hits: %d, Misses: %d, Size: %d, Efficiency: %.2f' % (self.__hits, self.__misses, len(self),
                                                                      self.__hits / self.__misses))

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

