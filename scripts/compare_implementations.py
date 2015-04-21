import os
import sys
import scipy.stats as stats

my_dir = os.path.dirname(os.path.abspath(__file__))

sys.path.append(os.path.join(my_dir, "../lib/python"))
from decaying_histogram import DecayingHistogram as c_DHist
from roll import *
from decay_equations import find_optimal_decay

if __name__ == '__main__':
    target_num_buckets = 50
    alpha_slow = 0.0001
    alpha_fast = find_optimal_decay(alpha_slow)
    c_slow = c_DHist(target_num_buckets, alpha_slow)
    c_fast = c_DHist(target_num_buckets, alpha_fast)

    py_slow = create_rolling_histogram_class(
            Bucket=create_bucket_class(alpha_count=alpha_slow),
            target_buckets=target_num_buckets)()
    py_fast = create_rolling_histogram_class(
            Bucket=create_bucket_class(alpha_count=alpha_fast),
            target_buckets=target_num_buckets)()

    dist = stats.norm(0, 1)
    num_iterations = 1000000
    for idx in range(num_iterations):
        if idx % int(num_iterations / 100) == 0:
            print >> sys.stderr, "{0} / {1}\r".format(idx, num_iterations),
            sys.stderr.flush()
        observation = dist.rvs(1)[0]
        #c_slow.insert(observation)
        #c_fast.insert(observation)
        py_slow.update(observation)
        py_fast.update(observation)
    for idx in range(num_iterations):
        if idx % int(num_iterations / 100) == 0:
            print >> sys.stderr, "{0} / {1}\r".format(idx, num_iterations),
            sys.stderr.flush()
        observation = dist.rvs(1)[0]
        #c_slow.insert(observation)
        #c_fast.insert(observation)
        #print c_slow.Jaccard_distance(c_fast)

        py_slow.update(observation)
        py_fast.update(observation)
        cdf_slow = py_slow.get_CDF()
        cdf_fast = py_fast.get_CDF()
        print calc_jaccard_distance(cdf_slow, cdf_fast)

    #cdf_slow = py_slow.get_CDF()
    #cdf_fast = py_fast.get_CDF()
    #print "====="
    #print "py: ", calc_jaccard_distance(cdf_slow, cdf_fast)
    #print "c:  ", c_slow.Jaccard_distance(c_fast)
    #print c_slow
    #print c_fast

