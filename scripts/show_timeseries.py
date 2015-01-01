from roll import *
import datasets

def simple_timeseries(data):
    for d in data:
        pylab.plot(range(len(d)), d, color="black", alpha=0.2)
    pylab.show()

if __name__ == '__main__':
    d = stats.norm(0, 1)

    ts_data = []
    size = 3000

    dists = [stats.expon(-2, 2), stats.norm(4, 2), stats.norm(-6, 1)]
    upper = 10
    for i in range(upper):
        print "Generating data... {0} / {1}\r".format(i, upper),
        sys.stdout.flush()
        data = []
        rh_slow = create_rolling_histogram_class(
                Bucket=create_bucket_class(
                    alpha_mu=0.005, alpha_count=0.005),
                target_buckets=200)()
        rh_fast = create_rolling_histogram_class(
                Bucket=create_bucket_class(
                    alpha_mu=0.005, alpha_count=0.01),
                target_buckets=75)()

        for val in datasets.dataset1():
            rh_slow.update(val)
            rh_fast.update(val)
            cdf_long = rh_slow.get_CDF()
            cdf_short = rh_fast.get_CDF()
            data.append(calc_jaccard_distance(cdf_long, cdf_short))

    #   #for val in gen_value(dists[:2], size):
    #   for val in gen_value([stats.norm(0, 1)], size):
    #       rh_slow.update(val)
    #       rh_fast.update(val)
    #       cdf_long = rh_slow.get_CDF()
    #       cdf_short = rh_fast.get_CDF()
    #       data.append(calc_KS(cdf_long, cdf_short))

    #   #for val in gen_value(dists[2:3], size):
    #   for val in gen_value([stats.norm(0, 2)], size):
    #       rh_slow.update(val)
    #       rh_fast.update(val)
    #       cdf_long = rh_slow.get_CDF()
    #       cdf_short = rh_fast.get_CDF()
    #       data.append(calc_KS(cdf_long, cdf_short))

    #   #for val in gen_value(dists[:2], size):
    #   for val in gen_value([stats.norm(0, 1)], size):
    #       rh_slow.update(val)
    #       rh_fast.update(val)
    #       cdf_long = rh_slow.get_CDF()
    #       cdf_short = rh_fast.get_CDF()
    #       data.append(calc_KS(cdf_long, cdf_short))

        ts_data.append(data)
    #simple_histogram(data)
    simple_timeseries(ts_data)

