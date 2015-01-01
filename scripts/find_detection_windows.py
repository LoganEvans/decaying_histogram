import multiprocessing
import pylab
from scipy import stats
import numpy as np
import sys
import os
import sys
from roll import create_rolling_histogram_class, create_bucket_class
from roll import gen_value, get_KS
import time

import decay_equations
import datasets
from show_density import ergodic_chain, simple_histogram


dist0 = stats.norm(0, 1)
#dist1 = stats.norm(0, 1.5)
dist1 = datasets.dataset1()


def get_detection_windows(args):
    (burnin, run_length, alpha_count_slow, alpha_count_fast,
     alpha_mu_slow, alpha_mu_fast, buckets_slow, buckets_fast,
     threshold) = args

    np.random.seed((os.getpid() << 16) | (int(time.time()) & 0xFFFF))
    rh_slow = create_rolling_histogram_class(
            Bucket=create_bucket_class(
                alpha_mu=alpha_mu_slow, alpha_count=alpha_count_slow),
            target_buckets=buckets_slow)()
    rh_fast = create_rolling_histogram_class(
            Bucket=create_bucket_class(
                alpha_mu=alpha_mu_fast, alpha_count=alpha_count_fast),
            target_buckets=buckets_fast)()

   #for val in gen_value([dist0], burnin):
   #    rh_slow.update(val)
   #    rh_fast.update(val)

    data = []
    #for i, val in enumerate(gen_value([dist1], run_length)):
    for i, val in enumerate(datasets.dataset1()):
        rh_slow.update(val)
        rh_fast.update(val)
        cdf_slow = rh_slow.get_CDF()
        cdf_fast = rh_fast.get_CDF()
        KS = get_KS(cdf_slow, cdf_fast)
        if KS > threshold:
            data.append(i)
    return data

def main():
    burnin = 10000
    significance_samples = 10000
    per_process_samples = significance_samples / multiprocessing.cpu_count()
    alpha_count_slow = 0.005
    alpha_count_fast = decay_equations.find_optimal_decay(alpha_count_slow)
    alpha_mu_slow = 0.01
    alpha_mu_fast = 0.01
    buckets_slow = 100
    buckets_fast = 100
    total_work = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    runs = pool.map(
            ergodic_chain,
            [[burnin, per_process_samples,
              alpha_count_slow, alpha_count_fast,
              alpha_mu_slow, alpha_mu_fast,
              buckets_slow, buckets_fast] for _ in range(total_work)])
    data = []
    for run in runs:
        for val in run:
            data.append(val)
    data.sort()
    threshold = data[int(0.9999 * len(data))]

    total_work = 100
    upper_bound = 10000
    windows = pool.map(
            get_detection_windows,
            [[burnin, upper_bound,
              alpha_count_slow, alpha_count_fast,
              alpha_mu_slow, alpha_mu_fast,
              buckets_slow, buckets_fast,
              threshold]
             for _ in range(total_work)])
    window_data = []
    for window in windows:
        for datum in window:
            window_data.append(datum)

    for _ in range(total_work):
        window_data.append(0)

    n, bins, patches = pylab.hist(
            window_data, upper_bound, normed=False, histtype='stepfilled',
            range=(0, upper_bound))
    pylab.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
    pylab.show()


if __name__ == '__main__':
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    main()

