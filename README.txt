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
$ src/profile_malloc -a -h 50 -t 3 -m 123KiB | scripts/animate.py
or
$ LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libjemalloc.so" src/profile_malloc -a -h 50 -t 3 | scripts/animate.py

The flag "-a" enables animation output, "-h 50" creates a decaying histogram
with a target of 50 buckets, the "-t 3" enables three threads, and "-m 123KiB"
instructs threads to malloc 123KiB at a time. See the code in
src/profile_malloc.cpp for other commandline options.

Tests
-----
To build the tests, run:
$ ./autogen.sh
$ ./configure
$ make check
$ tests/test_dhist

