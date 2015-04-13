import ctypes

class _DecayingHistogram(ctypes.Structure):
    _fields_ = [
            ("delete_bucket_threshold", ctypes.c_double),
            ("split_bucket_threshold", ctypes.c_double),
            ("alpha", ctypes.c_double),
            ("generation", ctypes.c_ulonglong),
            ("root", ctypes.c_void_p),
            ("bucket_list", ctypes.c_void_p),
            ("num_buckets", ctypes.c_uint),
            ("max_num_buckets", ctypes.c_uint),
            ("pow_table", ctypes.c_void_p),
            ("tree_mtx", ctypes.c_void_p),
            ("generation_mtx", ctypes.c_void_p)]

class DecayingHistogram(object):
    _lib = ctypes.CDLL("../../src/.libs/libdhistlib.so")
    def __init__(self, target_buckets, alpha):
        self._dhist = _DecayingHistogram()
        self._dhist_ptr = ctypes.pointer(self._dhist)
        self.DHIST_SINGLE_THREADED = ctypes.cast(
                self._lib.DHIST_SINGLE_THREADED,
                ctypes.POINTER(ctypes.c_int)).contents.value
        self.DHIST_MULTI_THREADED = ctypes.cast(
                self._lib.DHIST_MULTI_THREADED,
                ctypes.POINTER(ctypes.c_int)).contents.value
        self._lib.init_decaying_histogram(
                self._dhist_ptr, ctypes.c_uint(target_buckets),
                ctypes.c_double(alpha))

    def insert(self, observation):
        self._lib.dh_insert(
                self._dhist_ptr, ctypes.c_double(observation),
                self.DHIST_SINGLE_THREADED)

    def print_histogram(self):
        self._lib.print_histogram(
                self._dhist_ptr,
                ctypes.c_bool(True),
                ctypes.c_char_p("foobar"),
                ctypes.c_char_p("xaxis"),
                ctypes.c_char_p("yaxis"))


foo = DecayingHistogram(50, 0.0001)
import random
while True:
    for idx in range(10000):
        foo.insert(random.normalvariate(0, 1))
    foo.print_histogram()

