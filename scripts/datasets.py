import scipy.stats as stats

def dataset1():
    """Dataset 1 from
    http://sugiyama-www.cs.titech.ac.jp/~sugi/2009/SDM2009b.pdf

    """
    y = [0.0, 0.0]
    mu = 0.0
    for i in xrange(10):
        mu = mu + i
        noise = stats.norm(mu, 1)
        for _ in xrange(1000):
            y = [y[-1], 0.6 * y[-1] - 0.5 * y[-2] + noise.rvs(1)[0]]
            yield y[-1]


def dataset2():
    """Dataset 2 from
    http://sugiyama-www.cs.titech.ac.jp/~sugi/2009/SDM2009b.pdf

    """
    y = [0.0, 8.0, 6.0, 4.0]
    for val in y:
        yield val

    noise = stats.norm(0, 0.2)

    for _ in xrange(1000):
        # This is malarky.
        y = y[1:] + [noise.rvs(1)[0] + 0.97 * y[-1] + y[-2] -
                     0.5 * y[-3] + 0.97 * y[-4]]
        yield y[-1]

    for _ in xrange(500):
        y = y[1:] + [noise.rvs(1)[0] + 0.97 * y[-1] + y[-2] -
                     0.7 * y[-3] + 0.97 * y[-4]]
        yield y[-1]


if __name__ == '__main__':
    for v in dataset1():
        print v

