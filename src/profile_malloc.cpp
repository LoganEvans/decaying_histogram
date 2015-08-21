/*
 * Copyright (c) 2015, Logan P. Evans <loganpevans@gmail.com>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */


#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/types.h>

#include <atomic>
#include <random>

#include "dhist.h"

#define SUCCESS 9001

static struct dhist *g_histogram;
static std::atomic<char *> *g_buffer;
static uint64_t g_buffer_len;

static uint64_t rdtsc()
{
  uint32_t hi, lo;

  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return (((uint64_t)lo) | (((uint64_t)hi) << 32));
}

struct thread_func_args {
  uint64_t end_timestamp;
  uint64_t bytes_per_malloc;
  uint64_t num_mallocs_out;
  uint64_t total_cycles_out;
  bool profile_free;
  uint8_t __padding[7];
};

static void
bytes_to_display(uint64_t bytes, char *buffer) {
  if (bytes >= 1024 * 1024 * 1024) {
    sprintf(buffer, "%.3lfGiB", bytes / (1024 * 1024 * 1024.0));
  } else if (bytes >= 1024 * 1024) {
    sprintf(buffer, "%.3lfMiB", bytes / (1024 * 1024.0));
  } else if (bytes >= 1024) {
    sprintf(buffer, "%.3lfKiB", bytes / (1024.0));
  } else {
    sprintf(buffer, "%dB", bytes);
  }
}

static void *
thread_func(void *args) {
  uint64_t start_timestamp, stop_timestamp, end_timestamp;
  uint64_t bytes_per_malloc;
  uint64_t num_mallocs_out;
  uint64_t total_cycles_out;
  bool profile_free;
  int malloc_idx;
  std::uniform_int_distribution<int> *distribution;
  std::default_random_engine *generator;
  char *pointer;

  start_timestamp = stop_timestamp = 0;
  end_timestamp = ((struct thread_func_args *)args)->end_timestamp;
  bytes_per_malloc = ((struct thread_func_args *)args)->bytes_per_malloc;
  profile_free = ((struct thread_func_args *)args)->profile_free;
  num_mallocs_out = 0;
  total_cycles_out = 0;

  distribution =
      new std::uniform_int_distribution<int>(0, (int)(g_buffer_len - 1));
  generator = new std::default_random_engine((int)pthread_self());

  while (end_timestamp == 0 || rdtsc() < end_timestamp) {
    malloc_idx = (*distribution)(*generator);

    pointer = g_buffer[malloc_idx].exchange((char *)g_buffer);
    if (pointer == NULL) {
      // Nothing was here, so malloc and set.
      start_timestamp = rdtsc();
      g_buffer[malloc_idx].store((char *)malloc(bytes_per_malloc));
      stop_timestamp = rdtsc();
      if (g_histogram && !profile_free) {
        dhist_insert(
            g_histogram, log2(stop_timestamp - start_timestamp),
            DHIST_MULTI_THREADED);
      }
      num_mallocs_out++;
      total_cycles_out += stop_timestamp - start_timestamp;
    } else if (pointer == (char *)g_buffer) {
      // Some other thread was working on this item, so find some other index.
      continue;
    } else {
      // A pointer was here, so free it.
      if (profile_free)
        start_timestamp = rdtsc();
      free(pointer);
      if (profile_free)
        stop_timestamp = rdtsc();
      g_buffer[malloc_idx].store(NULL);

      if (profile_free) {
        dhist_insert(
            g_histogram, log2(stop_timestamp - start_timestamp),
            DHIST_MULTI_THREADED);
      }
    }
  }

  ((struct thread_func_args *)args)->num_mallocs_out = num_mallocs_out;
  ((struct thread_func_args *)args)->total_cycles_out = total_cycles_out;

  delete distribution;
  delete generator;
  return (void *)SUCCESS;
}

int main(int argc, char **argv) {
  pthread_t *threads;
  struct thread_func_args *args;
  struct timespec tim, tim2;
  char title[1000];
  char bytes_per_malloc_display[100];
  char num_bytes_display[100];
  char *histogram_json;
  uint64_t end_timestamp;
  uint64_t num_mallocs;
  uint64_t total_cycles;
  void *status;

  char ch;
  extern char *optarg;
  char optarg_copy[100];
  size_t len;
  double bytes;  // Units unknown.

  uint16_t opt_t_num_threads = 1;
  uint64_t opt_b_num_bytes = 1UL * 1024 * 1024 * 1024;
  uint32_t opt_m_bytes_per_malloc = 1024 * 1024;
  uint16_t opt_h_num_histogram_buckets = 0;
  double opt_d_decay_rate = 0.99999;
  double opt_s_seconds = 2.0;
  double opt_p_fps = 10.0;
  bool opt_f_profile_free = false;
  bool opt_a_animate = false;

  while ((ch = (char)getopt(argc, argv, "t:b:m:h:d:s:p:f:a")) != -1) {
    switch (ch) {
    case 't':
      // Number of threads.
      opt_t_num_threads = (uint16_t)strtoul(optarg, NULL, 10);
      break;
    case 'b':
      // Permitted number of bytes, shared among all threads.
      opt_b_num_bytes = strtoull(optarg, NULL, 10);
      break;
    case 'm':
      // malloc size in bytes.

      len = strlen(optarg);
      if (len <= strlen("..")) {
        opt_m_bytes_per_malloc = (uint32_t)strtoul(optarg, NULL, 10);
      } else {
        strncpy(optarg_copy, optarg, sizeof(optarg_copy));
        if (strcmp(optarg + len - strlen("GiB"), "GiB") == 0) {
          optarg_copy[len - strlen("GiB")] = '\0';
          bytes = strtod(optarg_copy, NULL);
          opt_m_bytes_per_malloc = (uint32_t)(bytes * 1024 * 1024 * 1024UL);
        } else if (strcmp(optarg + len - strlen("MiB"), "MiB") == 0) {
          optarg_copy[len - strlen("MiB")] = '\0';
          bytes = strtod(optarg_copy, NULL);
          opt_m_bytes_per_malloc = (uint32_t)(bytes * 1024 * 1024UL);
        } else if (strcmp(optarg + len - strlen("KiB"), "KiB") == 0) {
          optarg_copy[len - strlen("KiB")] = '\0';
          bytes = strtod(optarg_copy, NULL);
          opt_m_bytes_per_malloc = (uint32_t)(bytes * 1024UL);
        } else {
          opt_m_bytes_per_malloc = (uint32_t)strtoul(optarg, NULL, 10);
        }
      }
      break;
    case 'h':
      // If present, use a decaying histogram with the specified number of
      // buckets.
      opt_h_num_histogram_buckets = (uint16_t)strtoul(optarg, NULL, 10);
      break;
    case 'd':
      // Decay rate for the decaying histogram.
      opt_d_decay_rate = strtod(optarg, NULL);
      break;
    case 's':
      // Rough number of seconds (based on a 2.3 GHz cpu).
      opt_s_seconds = strtod(optarg, NULL);
      break;
    case 'p':
      // FPS.
      opt_p_fps = strtod(optarg, NULL);
      break;
    case 'f':
      // Profile free instead of malloc.
      opt_f_profile_free = true;
      break;
    case 'a':
      // Animate. If this is present, -s is ignored.
      opt_a_animate = true;
      break;
    default:
      break;
    }
  }

  if (!opt_h_num_histogram_buckets && (opt_a_animate || opt_f_profile_free)) {
    fprintf(stderr, "ERROR: -a and -f options depend on -h\n");
    assert(false);
  }

  if (opt_h_num_histogram_buckets)
    g_histogram = dhist_init(opt_h_num_histogram_buckets, opt_d_decay_rate);
  else
    g_histogram = NULL;

  if (!g_histogram)
      printf("threads,cycles_per_malloc,num_mallocs,bytes_per_malloc\n");

  if (opt_h_num_histogram_buckets) {
    bytes_to_display(opt_m_bytes_per_malloc, bytes_per_malloc_display);
    bytes_to_display(opt_b_num_bytes, num_bytes_display);

    sprintf(title,
        "threads: %d, size-per-malloc: %s, bytes-permitted: %s",
        opt_t_num_threads, bytes_per_malloc_display, num_bytes_display);
  }

  g_buffer_len = opt_b_num_bytes / opt_m_bytes_per_malloc;
  g_buffer = (std::atomic<char *>*)malloc(
      sizeof(std::atomic<char *>) * g_buffer_len);
  for (uint64_t idx = 0; idx < g_buffer_len; idx++)
    g_buffer[idx].store(NULL);

  threads = (pthread_t *)malloc(sizeof(pthread_t) * opt_t_num_threads);
  args = (struct thread_func_args *)malloc(
      sizeof(struct thread_func_args) * opt_t_num_threads);
  if (opt_a_animate) {
    end_timestamp = 0;
  } else {
    end_timestamp =
        rdtsc() + (uint64_t)(opt_s_seconds * 2.3 * 1024 * 1024 * 1024);
  }

  for (int i = 0; i < opt_t_num_threads; i++) {
    args[i].end_timestamp = end_timestamp;
    args[i].bytes_per_malloc = opt_m_bytes_per_malloc;
    args[i].profile_free = opt_f_profile_free;
    pthread_create(
        &threads[i], NULL, (void *(*)(void *))thread_func, &args[i]);
  }

  if (opt_a_animate) {
    tim.tv_sec = (int)(1.0 / opt_p_fps);
    tim.tv_nsec = (int)(1000000000 * (1.0 / opt_p_fps - tim.tv_sec));

    while (1) {
      nanosleep(&tim, &tim2);
      histogram_json = dhist_get_json(
          g_histogram, title, "log_2(insertion time)",
          DHIST_MULTI_THREADED);
      puts(histogram_json);
      free(histogram_json);
      fflush(stdout);
    }
  }

  num_mallocs = 0;
  total_cycles = 0;
  for (int i = 0; i < opt_t_num_threads; i++) {
    pthread_join(threads[i], &status);
    if ((uint64_t)status != SUCCESS) {
      printf("Thread %d returned status: %lu\n", i, (uint64_t)status);
      assert(false);
    }
    num_mallocs += args[i].num_mallocs_out;
    total_cycles += args[i].total_cycles_out;
  }

  free(threads);
  free(args);
  for (uint64_t idx = 0; idx < g_buffer_len; idx++)
    free(g_buffer[idx].load());
  free(g_buffer);

  if (g_histogram) {
    if (!opt_a_animate) {
      histogram_json = dhist_get_json(
          g_histogram, title, "log_2(insertion time)",
          DHIST_MULTI_THREADED);
      puts(histogram_json);
      free(histogram_json);
    }
    dhist_destroy(g_histogram);
  } else {
    printf("%d,%lu,%lu,%u,%lu\n",
        opt_t_num_threads, total_cycles / num_mallocs, num_mallocs,
        opt_m_bytes_per_malloc, opt_b_num_bytes / opt_t_num_threads);
  }

  return 0;
}

