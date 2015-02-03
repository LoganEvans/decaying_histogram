To build:
$ autoreconf -ivf
$ ./configure CC=clang CXX=clang++ CFLAGS="-O2 -g" CXXFLAGS="-O2 -g"
$ make

To build the tests, first ./configure, and then:
$ make check TESTS=
$ tests/decaying_histogram_test

To see the decaying histogram at work, try:
$ src/main | scripts/animate.py

