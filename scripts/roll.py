from Queue import Queue
import numpy as np
import scipy.stats as stats
import sys
import time
from pprint import pprint

def create_bucket_class(alpha_mu=0.001, alpha_count=0.001):
    class Bucket(object):
        ALPHA_MU = None
        ALPHA_COUNT = None

        def __init__(self, mu, count, histogram):
            self.mu = mu
            self.count = count
            self.bucket_below = None
            self.bucket_above = None
            self._lower_bound = None
            self._upper_bound = None
            # This is to avoid a ZeroDivisionError
            self._max = mu + 0.0001
            self._min = mu
            self._last_decay_generation = 0
            self.histogram = histogram

        def is_in_bucket(self, observation):
            if ((self._lower_bound is None or self._lower_bound <= observation) and
                (self._upper_bound is None or observation <= self._upper_bound)
            ):
                return True
            else:
                return False

        @staticmethod
        def recompute_bound(lower, upper):
            """Both buckets should be decayed."""
            assert lower.bucket_above == upper
            assert upper.bucket_below == lower
            bound = ((lower.mu * lower.count + upper.mu * upper.count) /
                     (lower.count + upper.count))
            lower._upper_bound = bound
            upper._lower_bound = bound

        def update(self, observation):
            assert self.is_in_bucket(observation)
            if observation > self._max:
                self._max = observation
            if observation < self._min:
                self._min = observation

            self.decay()
            self.count += 1.0
            self.mu = (
                    Bucket.ALPHA_MU * observation +
                    (1.0 - Bucket.ALPHA_MU) * self.mu)

            if self.bucket_below:
                self.bucket_below.decay()
                Bucket.recompute_bound(self.bucket_below, self)

            if self.bucket_above:
                self.bucket_above.decay()
                Bucket.recompute_bound(self, self.bucket_above)

        def lower_bound(self):
            if self._lower_bound is None:
                return self._min
            else:
                return self._lower_bound

        def upper_bound(self):
            if self._upper_bound is None:
                return self._max
            else:
                return self._upper_bound

        def decay(self):
            decay_generation = self.histogram.get_decay_generation()
            if decay_generation == self._last_decay_generation:
                return
            self.count = self.count * (
                    (1.0 - Bucket.ALPHA_COUNT) **
                    (decay_generation - self._last_decay_generation))
            self._last_decay_generation = decay_generation

        def density(self):
            self.decay()
            return Bucket.ALPHA_COUNT * self.count

        def get_height(self):
            return self.density() / (self.upper_bound() - self.lower_bound())

    Bucket.ALPHA_MU = alpha_mu
    Bucket.ALPHA_COUNT = alpha_count
    return Bucket

def create_rolling_histogram_class(Bucket=create_bucket_class(),
                                   target_buckets=100):
    class RollingHistogram(object):
        TARGET_BUCKETS = None

        def __init__(self):
            self.expected_count = (
                    1.0 /
                    (Bucket.ALPHA_COUNT * RollingHistogram.TARGET_BUCKETS))
            self.delete_bucket_threshold = self.expected_count / 2.0
            self.split_bucket_threshold = self.expected_count * 2.0

            # While we can start with one bucket, the histogram code will
            # complain if we have fewer than two.
            self.bucket_list = [
                    Bucket(0.0, self.expected_count, self),
                    Bucket(1.0, self.expected_count, self),
                ]
            self.bucket_list[0].bucket_above = self.bucket_list[1]
            self.bucket_list[1].bucket_below = self.bucket_list[0]
            self._decay_generation = 0

        def get_decay_generation(self):
            return self._decay_generation

        def find_bucket(self, observation):
            lo = 0
            hi = len(self.bucket_list)
            while lo < hi:
                mid = (lo + hi) / 2
                midval = self.bucket_list[mid]
                if self.bucket_list[mid].is_in_bucket(observation):
                    return mid
                elif self.bucket_list[mid].mu < observation:
                    lo = mid + 1
                else:
                    hi = mid
            assert False

        def full_refresh(self):
            for bucket in self.bucket_list:
                bucket.decay()

            actions = [
                    [lambda count: count < self.delete_bucket_threshold,
                     self.delete_bucket],
                    [lambda count: self.split_bucket_threshold < count,
                     self.split_bucket]]
            do_another_round = True
            # We want to do all of the splits and then all of the deletes
            while do_another_round:
                do_another_round = False
                for action_pair in actions:
                    test, action = action_pair
                    bucket_idx = 0
                    while bucket_idx < len(self.bucket_list):
                        count = self.bucket_list[bucket_idx].count
                        if test(count):
                            action(bucket_idx)
                            do_another_round = True
                            continue
                        bucket_idx += 1

        def update(self, observation):
            self._decay_generation += 1

            bucket_idx = self.find_bucket(observation)
            bucket = self.bucket_list[bucket_idx]
            bucket.update(observation)

            if bucket.count < self.delete_bucket_threshold:
                self.full_refresh()
            if bucket.count > self.split_bucket_threshold:
                self.full_refresh()

        def get_histogram(self):
            n = [bucket.get_height() for bucket in self.bucket_list]
            bins = [bucket.upper_bound() for bucket in self.bucket_list]
            bins[-1] = (2 * self.bucket_list[-1].mu -
                        self.bucket_list[-1].lower_bound())
            return np.array(n), np.array(bins)

        def get_CDF(self):
            self.full_refresh()
            cdf = []
            acc = 0.0
            for bucket in self.bucket_list:
                acc += bucket.density()
                cdf.append([bucket.upper_bound(), acc])
            return cdf

        def density(self, observation):
            bucket_idx = self.find_bucket(observation)
            return self.bucket_list[bucket_idx].density()

        def delete_bucket(self, bucket_idx):
            self.bucket_list[bucket_idx].decay()
            count_below = None
            count_above = None
            if bucket_idx - 1 >= 0:
                self.bucket_list[bucket_idx - 1].bucket_above = (
                        self.bucket_list[bucket_idx].bucket_above)
                self.bucket_list[bucket_idx - 1].decay()
                count_below = self.bucket_list[bucket_idx - 1].count

            if bucket_idx + 1 < len(self.bucket_list):
                self.bucket_list[bucket_idx + 1].bucket_below = (
                        self.bucket_list[bucket_idx].bucket_below)
                self.bucket_list[bucket_idx + 1].decay()
                count_above = self.bucket_list[bucket_idx + 1].count

            if None not in [count_below, count_above]:
                if count_below < count_above:
                    lucky_bucket = self.bucket_list[bucket_idx - 1]
                else:
                    lucky_bucket = self.bucket_list[bucket_idx + 1]
            elif count_below is not None:
                lucky_bucket = self.bucket_list[bucket_idx - 1]
            else:
                assert len(self.bucket_list) > 1
                lucky_bucket = self.bucket_list[bucket_idx + 1]

            dying_bucket = self.bucket_list[bucket_idx]
            lucky_bucket.mu = ((
                    lucky_bucket.mu * lucky_bucket.count +
                    dying_bucket.mu * dying_bucket.count) /
                    (lucky_bucket.count + dying_bucket.count))
            lucky_bucket.count += dying_bucket.count
            lucky_bucket._max = max(lucky_bucket._max, dying_bucket._max)
            lucky_bucket._min = min(lucky_bucket._min, dying_bucket._min)
            self.bucket_list.pop(bucket_idx)
            if lucky_bucket.bucket_below:
                lucky_bucket.bucket_below.decay()
                Bucket.recompute_bound(lucky_bucket.bucket_below, lucky_bucket)
            else:
                lucky_bucket._lower_bound = None

            if lucky_bucket.bucket_above:
                lucky_bucket.bucket_above.decay()
                Bucket.recompute_bound(lucky_bucket, lucky_bucket.bucket_above)
            else:
                lucky_bucket._upper_bound = None

        def assert_things(self, idx):
            show = False
            for bucket in self.bucket_list[1:-1]:
                if bucket.bucket_below._upper_bound > bucket._upper_bound:
                    show = True
                if bucket._lower_bound > bucket.bucket_above._lower_bound:
                    show = True
            if show:
                print
                print idx
                for bucket in self.bucket_list:
                    print "{0:>8.4} | {1:>8.4} {2}".format(
                            bucket._lower_bound, bucket._upper_bound,
                            bucket.count)
                assert False

        def split_bucket(self, bucket_idx):
            bucket = self.bucket_list[bucket_idx]
            lower = bucket.lower_bound()
            upper = bucket.upper_bound()
            diameter = upper - lower
            median = lower + diameter / 2.0
            mu_upper = median + diameter / 6.0
            mu_lower = median - diameter / 6.0
            count_upper = bucket.count / 2.0
            count_lower = bucket.count / 2.0

            new_bucket = Bucket(mu_lower, count_lower, self)
            new_bucket._last_decay_generation = bucket._last_decay_generation
            bucket.mu = mu_upper
            bucket.count = count_upper

            if bucket.bucket_below:
                bucket.bucket_below.bucket_above = new_bucket
                new_bucket.bucket_below = bucket.bucket_below
                Bucket.recompute_bound(new_bucket.bucket_below, new_bucket)

            new_bucket.bucket_above = bucket
            bucket.bucket_below = new_bucket
            Bucket.recompute_bound(new_bucket, new_bucket.bucket_above)

            if bucket.bucket_above:
                Bucket.recompute_bound(bucket, bucket.bucket_above)

            self.bucket_list.insert(bucket_idx, new_bucket)

        def display(self):
            from matplotlib import pyplot
            bins = [bucket.lower_bound() for bucket in self.bucket_list]
            bins.append(self.bucket_list[-1].upper_bound())
            bins = np.array(bins)
            widths = bins[1:] - bins[:-1]
            freqs = np.array([bucket.density() for bucket in self.bucket_list])
            heights = freqs.astype(np.float) / widths
            pyplot.fill_between(
                    bins.repeat(2)[1:-1], heights.repeat(2),
                    facecolor='steelblue')
            pyplot.show()

    RollingHistogram.TARGET_BUCKETS = target_buckets
    return RollingHistogram


def redistribute(cdf, new_cdf):
    insert_idx = 0
    for i in range(len(cdf) - 1):
        slope = (cdf[i + 1][1] - cdf[i][1]) / (cdf[i + 1][0] - cdf[i][0])
        while new_cdf[insert_idx][0] < cdf[i + 1][0]:
            if insert_idx == len(new_cdf):
                break
            new_cdf[insert_idx][1] = (
                    cdf[i][1] + slope * (new_cdf[insert_idx][0] - cdf[i][0]))
            insert_idx += 1
    new_cdf[-1][1] = 1.0


def calc_jaccard_distance(cdf0, cdf1):
    x_values = sorted(set(zip(*cdf0)[0] + zip(*cdf1)[0]))
    new_cdf0 = [[x, None] for x in x_values]
    new_cdf1 = [[x, None] for x in x_values]

    redistribute(cdf0, new_cdf0)
    redistribute(cdf1, new_cdf1)

    intersection = 0.0
    union = 0.0
    for i in range(len(new_cdf0)):
        intersection += min(new_cdf0[i][1], new_cdf1[i][1])
        union += max(new_cdf0[i][1], new_cdf1[i][1])
    return 1.0 - intersection / union


def calc_cramer_von_mises_criterion(cdf0, cdf1):
    """XXX This seems wrong..."""
    x_values = sorted(set(zip(*cdf0)[0] + zip(*cdf1)[0]))
    new_cdf0 = [[x, None] for x in x_values]
    new_cdf1 = [[x, None] for x in x_values]

    redistribute(cdf0, new_cdf0)
    redistribute(cdf1, new_cdf1)

    g = lambda t: np.abs(new_cdf0[t][1] - new_cdf1[t][1]) ** 3 / 3.0
    acc = 0.0
    for i in range(len(new_cdf0)):
        acc += g(i)
    return acc


def calc_kullback_leibler_divergence(cdf0, cdf1):
    x_values = sorted(set(zip(*cdf0)[0] + zip(*cdf1)[0]))
    new_cdf0 = [[x, None] for x in x_values]
    new_cdf1 = [[x, None] for x in x_values]
    dist0 = [None for x in range(len(x_values) - 1)]
    dist1 = [None for x in range(len(x_values) - 1)]

    redistribute(cdf0, new_cdf0)
    redistribute(cdf1, new_cdf1)

    for i in range(1, len(new_cdf0)):
        dist0[i - 1] = new_cdf0[i][1] - new_cdf0[i - 1][1]
        dist1[i - 1] = new_cdf1[i][1] - new_cdf1[i - 1][1]

    return stats.entropy(dist0, dist1)

def calc_KS(cdf0, cdf1):
    idx0, idx1 = 0, 0
    KS = 0.0
    partial0 = 0.0
    partial1 = 0.0

    len0 = len(cdf0)
    len1 = len(cdf1)
    while idx0 < len0 and idx1 < len1:
        KS = max(KS, np.abs(cdf0[idx0][1] - cdf1[idx1][1]))
        if idx0 + 2 == len0 or idx1 + 2 == len1:
            # Since we don't record the right-hand values for the last bucket,
            # we can't comput that distance.
            break
        elif cdf0[idx0 + 1][0] == cdf1[idx1 + 1][0]:
            # The x values are exactly equal.
            idx0 += 1
            partial0 = 0.0
            idx1 += 1
            partial1 = 0.0
        elif cdf0[idx0 + 1][0] < cdf1[idx1 + 1][0]:
            # The x for 0 is smaller than for 1, so move the x value for 0.
            idx0 += 1
            partial0 = 0.0
            fraction = ((cdf0[idx0][0] - cdf1[idx1][0]) /
                        (cdf1[idx1 + 1][0] - cdf1[idx1][0]))
            partial1 = (cdf1[idx1 + 1][1] - cdf1[idx1][1]) * fraction
        else:
            # The x value for 1 is smaller than for 0, so move the x value for
            # 1.
            idx1 += 1
            partial1 = 0.0
            try:
                fraction = ((cdf1[idx1][0] - cdf0[idx0][0]) /
                            (cdf0[idx0 + 1][0] - cdf0[idx0][0]))
            except:
                from pprint import pprint
                pprint(cdf0)
                pprint(cdf1)
                raise
            partial0 = (cdf0[idx0 + 1][1] - cdf0[idx0][1]) * fraction
    return KS


def gen_value(dists, size):
    msize = 128
    selector_dist = stats.randint(0, len(dists))
    selector_Q = Queue(msize)
    dists_Q = [Queue(msize) for _ in range(len(dists))]
    for _ in xrange(size):
        if selector_Q.empty():
            to_add = selector_dist.rvs(msize)
            for val in to_add:
                selector_Q.put(val)

        selection = selector_Q.get()
        dist_Q = dists_Q[selection]
        if dist_Q.empty():
            to_add = dists[selection].rvs(msize)
            for val in to_add:
                dist_Q.put(val)

        yield dist_Q.get()


if __name__ == '__main__':
    d = stats.norm(0, 1)

    ts_data = []
    size = 5000

    dists = [stats.expon(-2, 2), stats.norm(4, 2), stats.norm(-6, 1)]
    for i in range(1):
        print i, '\r',
        sys.stdout.flush()
        data = []
        rh_slow = create_rolling_histogram_class(
                Bucket=create_bucket_class(
                    alpha_mu=0.005, alpha_count=0.001),
                target_buckets=50)()
        rh_fast = create_rolling_histogram_class(
                Bucket=create_bucket_class(
                    alpha_mu=0.005, alpha_count=0.005),
                target_buckets=50)()

        for val in gen_value(dists[:2], size):
            rh_slow.update(val)
            rh_fast.update(val)

        cdf_long = rh_slow.get_CDF()
        cdf_short = rh_fast.get_CDF()
        #print calc_jaccard_distance(cdf_short, cdf_long)
        #print calc_cramer_von_mises_criterion(cdf_short, cdf_long)
        print calc_kullback_leibler_divergence(cdf_short, cdf_long)

    #       cdf_long = rh_slow.get_CDF()
    #       cdf_short = rh_fast.get_CDF()
    #       #data.append(calc_KS(cdf_long, cdf_short))

    #   for val in gen_value(dists[2:], size):
    #       rh_slow.update(val)
    #       rh_fast.update(val)
    #       cdf_long = rh_slow.get_CDF()
    #       cdf_short = rh_fast.get_CDF()
    #       #data.append(calc_KS(cdf_long, cdf_short))

    #   for val in gen_value(dists[:2], size):
    #       rh_slow.update(val)
    #       rh_fast.update(val)
    #       cdf_long = rh_slow.get_CDF()
    #       cdf_short = rh_fast.get_CDF()
    #       #data.append(calc_KS(cdf_long, cdf_short))

        #ts_data.append(data)
    #simple_histogram(data)
    #simple_timeseries(ts_data)

