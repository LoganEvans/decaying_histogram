from matplotlib import pyplot
import numpy as np
import sys

freqs, bins = sys.argv[1].splitlines()
freqs = np.array(eval(freqs))
bins = np.array(eval(bins))
widths = bins[1:] - bins[:-1]
heights = freqs.astype(np.float)/widths
pyplot.fill_between(
        bins.repeat(2)[1:-1], heights.repeat(2), facecolor='steelblue')
pyplot.show()

