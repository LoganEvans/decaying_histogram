from matplotlib import pyplot
import numpy as np
import sys
import json

print '===='
print sys.argv[1]
histograms = [json.loads(line.strip()) for line in sys.argv[1].splitlines()]
for histogram in histograms:
    weights = np.array(histogram["weights"])
    bins = np.array(histogram["boundaries"])
    widths = bins[1:] - bins[:-1]
    heights = weights.astype(np.float) / widths
    pyplot.fill_between(
            bins.repeat(2)[1:-1], heights.repeat(2),
            color="steelblue", alpha=1.0 / len(histograms))
pyplot.show()

