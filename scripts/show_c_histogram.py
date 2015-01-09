from matplotlib import pyplot
import numpy as np
import sys
import json

data = json.loads(sys.argv[1])
densities = np.array(data['densities'])
boundaries = np.array(data['boundaries'])
widths = boundaries[1:] - boundaries[:-1]
heights = densities.astype(np.float)/widths
pyplot.fill_between(
        boundaries.repeat(2)[1:-1], heights.repeat(2), facecolor='steelblue')
pyplot.show()

