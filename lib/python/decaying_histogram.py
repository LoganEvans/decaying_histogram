import ctypes
import ctypes.util
import os
import json

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
    lib_dir = os.path.dirname(os.path.abspath(__file__))
    _libdhist = ctypes.CDLL(os.path.join(lib_dir, "libdhistlib.so"))
    _libdhist.get_new_histogram_json.argtypes = (
            ctypes.POINTER(_DecayingHistogram),
            ctypes.c_bool,
            ctypes.c_char_p,
            ctypes.c_char_p)
    # Not setting this to c_char_p because we have to free it.
    _libdhist.get_new_histogram_json.restype = ctypes.c_void_p

    _libc = ctypes.CDLL(ctypes.util.find_library('c'))
    _libc.free.argtypes = (ctypes.c_void_p,)

    def __init__(self, target_buckets, alpha):
        self._dhist = _DecayingHistogram()
        self._dhist_ptr = ctypes.pointer(self._dhist)
        self.DHIST_SINGLE_THREADED = ctypes.cast(
                self._libdhist.DHIST_SINGLE_THREADED,
                ctypes.POINTER(ctypes.c_int)).contents.value
        self.DHIST_MULTI_THREADED = ctypes.cast(
                self._libdhist.DHIST_MULTI_THREADED,
                ctypes.POINTER(ctypes.c_int)).contents.value
        self._libdhist.init_decaying_histogram(
                self._dhist_ptr, ctypes.c_uint(target_buckets),
                ctypes.c_double(alpha))

    def insert(self, observation):
        self._libdhist.dh_insert(
                self._dhist_ptr, ctypes.c_double(observation),
                self.DHIST_MULTI_THREADED)

    def get_json(self):
        return json.loads(str(self))

    def __str__(self):
        result = self._libdhist.get_new_histogram_json(
                self._dhist, True, "foobar", "xaxix", "yaxis")
        histogram_json = ctypes.cast(result, ctypes.c_char_p).value
        self._libc.free(result)
        return histogram_json


foo = DecayingHistogram(10, 0.0001)
import random
while True:
    for idx in range(10000):
        foo.insert(random.normalvariate(0, 1))
    x = foo.get_json()
    x["title"] = str(x["generation"])
    print json.dumps(x)

