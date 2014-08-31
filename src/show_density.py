import multiprocessing
import pylab
from scipy import stats
import numpy as np
import sys
import os
from roll import create_rolling_histogram_class, create_bucket_class
from roll import gen_value
from roll import calc_KS, calc_jaccard_distance
from roll import calc_cramer_von_mises_criterion
from roll import calc_kullback_leibler_divergence
from decay_equations import find_optimal_decay
import time

FUNC_LIST = [
        calc_KS, calc_jaccard_distance,
        calc_kullback_leibler_divergence,
        lambda cdf0, cdf1: calc_kullback_leibler_divergence(cdf1, cdf0)]
FUNC_LABELS = ["KS", "Jaccard", "KL short -> long", "KL long -> short"]

def simple_histogram(data):
    # the histogram of the data with histtype='step'
    n, bins, patches = pylab.hist(data, 250, normed=True, histtype='stepfilled')
    pylab.setp(patches, 'facecolor', 'g', 'alpha', 0.75)
    pylab.show()

def ergodic_chain(args):
    (burnin, run_length, alpha_count_slow, alpha_count_fast,
     alpha_mu_slow, alpha_mu_fast, buckets_slow, buckets_fast) = args

    np.random.seed((os.getpid() << 16) | (int(time.time()) & 0xFFFF))
    rh_slow = create_rolling_histogram_class(
            Bucket=create_bucket_class(
                alpha_mu=alpha_mu_slow, alpha_count=alpha_count_slow),
            target_buckets=buckets_slow)()
    rh_fast = create_rolling_histogram_class(
            Bucket=create_bucket_class(
                alpha_mu=alpha_mu_fast, alpha_count=alpha_count_fast),
            target_buckets=buckets_fast)()

    jagged = [stats.uniform(x, x + 1) for x in range(200)]
    #for val in gen_value(jagged, burnin):
    #for val in gen_value([stats.uniform(0, 1)], burnin):
    for val in gen_value([stats.norm(0, 1)], burnin):
    #for val in gen_value([stats.cauchy(0)], burnin):
        rh_slow.update(val)
        rh_fast.update(val)
    data = [[] for _ in range(len(FUNC_LIST))]
    #for val in gen_value(jagged, run_length):
    #for val in gen_value([stats.uniform(0, 1)], run_length):
    for val in gen_value([stats.norm(0, 1)], run_length):
    #for val in gen_value([stats.cauchy(0)], run_length):
        rh_slow.update(val)
        rh_fast.update(val)
        cdf_long = rh_slow.get_CDF()
        cdf_short = rh_fast.get_CDF()
        for i, func in enumerate(FUNC_LIST):
            data[i].append(func(cdf_short, cdf_long))
    return data


def main():
    total_work = multiprocessing.cpu_count()
    burnin = 30000
    significance_samples = 1000000
    per_process_samples = significance_samples / multiprocessing.cpu_count()
    alpha_count_slow = 0.001
    alpha_count_fast = find_optimal_decay(alpha_count_slow)
    alpha_mu_slow = 0.01
    alpha_mu_fast = 0.01
    buckets_slow = 50
    buckets_fast = 50
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    runs = pool.map(
            ergodic_chain,
            [[burnin, per_process_samples,
              alpha_count_slow, alpha_count_fast,
              alpha_mu_slow, alpha_mu_fast,
              buckets_slow, buckets_fast] for _ in range(total_work)])
    aggregator = [[] for _ in range(len(FUNC_LIST))]
    for run in runs:
        for i, data_list in enumerate(run):
            aggregator[i] += data_list
    colors = ['red', 'green', 'blue', 'purple']
    for label, data in zip(FUNC_LABELS, aggregator):
        #data.sort()
        _, _, patches = pylab.hist(
                data, 250, label=label,
                normed=True, histtype='stepfilled')
        pylab.setp(patches, 'alpha', 0.4)
    pylab.legend()
    pylab.show()
    #print np.mean(data), len(data)
    #print data[int(0.5 * len(data))], int(0.5 * len(data))
    #print data[int(0.9 * len(data))], int(0.9 * len(data))
    #print data[int(0.99 * len(data))], int(0.99 * len(data))
    #print data[int(0.999 * len(data))], int(0.999 * len(data))
    #print data[-1]

    #simple_histogram(data)

if __name__ == '__main__':
    #chart()
    main()

