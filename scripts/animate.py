#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot
import matplotlib.animation as animation
import sys
import json


_FRAMES = 10

def update(num):
    new_data = json.loads(sys.stdin.readline())
    update.history.insert(0, new_data)
    if len(update.history) > _FRAMES:
        update.history.pop()

    pyplot.cla()
    if "title" in new_data:
        pyplot.title(new_data["title"])
    if "xlabel" in new_data:
        pyplot.xlabel(new_data["xlabel"])
    for idx, data in enumerate(reversed(update.history)):
        weights = np.array(data['weights'])
        boundaries = np.array(data['boundaries'])
        widths = boundaries[1:] - boundaries[:-1]
        heights = weights.astype(np.float)/widths
        pyplot.fill_between(
                boundaries.repeat(2)[1:-1], heights.repeat(2),
                facecolor='steelblue',
                alpha=(idx + 1.0) / len(update.history))

update.history = []

fig1 = pyplot.figure()

pyplot.xlabel('x')
pyplot.ylabel('density')
pyplot.title('test')
line_ani = animation.FuncAnimation(fig1, update)

pyplot.show()


