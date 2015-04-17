AM_CXXFLAGS = $(WARN_CXXFLAGS)
AM_LDFLAGS = -pthread -lm

bin_PROGRAMS = main profile_malloc
lib_LTLIBRARIES = libdhistlib.la

libdhistlib_la_SOURCES = decaying_histogram.c

main_SOURCES = main.cpp
main_LDADD = libdhistlib.la
main_LDFLAGS = -no-install

profile_malloc_SOURCES = profile_malloc.cpp
profile_malloc_LDADD = libdhistlib.la
profile_malloc_LDFLAGS = -no-install


