To build
--------
$ ./autogen.sh
$ ./configure
$ make


Example
-------
To see the decaying histogram at work, try:
$ src/main | scripts/animate.py
This measures how much time it takes to insert an observation into a decaying
histogram and then inserts that latency into the histogram. In effect, this
profiles the decaying histogram code itself. By default, this will run with
4 threads (see the #defines in src/main.cpp).

Another example would be to explore the scaling characteristics of malloc:
$ src/profile_malloc | scripts/animate.py
or
$ LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libjemalloc.so" src/profile_malloc | scripts/animate.py


Tests
-----
To build the tests, first ./configure, and then:
$ make check TESTS=
$ tests/decaying_histogram_test

