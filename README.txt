To build:
$ autoreconf -ivf
$ ./configure CC=clang CXX=clang++ CFLAGS="-O2 -g" CXXFLAGS="-O2 -g"
$ make

To build the tests, first ./configure, and then:
$ make check TESTS=
$ tests/decaying_histogram_test

To see the decaying histogram at work, try:
$ src/main | scripts/animate.py

Another example would be to explore the scaling characteristics of malloc:
$ src/profile_malloc | scripts/animate.py
or
$ LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libjemalloc.so" src/profile_malloc | scripts/animate.py



