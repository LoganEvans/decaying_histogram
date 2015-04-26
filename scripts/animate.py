#!/usr/bin/env python

import argparse
import json
import numpy as np
import sys

from matplotlib import pyplot
import matplotlib.animation as animation

_color_order = [
        "blue", "green", "red", "cyan", "magenta", "yellow", "black"]

def parse_options():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '--frames', type=int, default=10)
    parser.add_argument(
            '--alpha', type=float, default=1.0,
            help="alpha transparency value for front histogram")

    return parser.parse_args()

def update(num, cli_args):
    global _color_order
    frames = cli_args.frames
    new_data = json.loads(sys.stdin.readline())
    if new_data["id"] not in update.memos:
        update.memos[new_data["id"]] = {
                "color": _color_order[0],
                "history": []}
        _color_order = _color_order[1:] + [_color_order[0]]
    new_memo = update.memos[new_data["id"]]

    if len(new_memo["history"]) > frames:
        new_memo["history"].pop()
    new_memo["history"].insert(0, new_data)

    pyplot.cla()
    if "title" in new_data:
        pyplot.title(new_data["title"])
    if "xlabel" in new_data:
        pyplot.xlabel(new_data["xlabel"])

    ymax = None
    for memo in update.memos.values():
        for idx, data in enumerate(reversed(memo["history"])):
            weights = np.array(data['weights'])
            boundaries = np.array(data['boundaries'])
            widths = boundaries[1:] - boundaries[:-1]
            max_weights = max(weights)
            if ymax is None or max_weights > ymax:
                ymax = max_weights
            pyplot.fill_between(
                    boundaries.repeat(2)[1:-1], weights.repeat(2),
                    facecolor=memo["color"],
                    alpha=cli_args.alpha * ((idx + 1.0) / len(memo["history"])))
    pyplot.ylim(0, ymax)
update.memos = {}

if __name__ == '__main__':
    cli_args = parse_options()

    fig1 = pyplot.figure()

    pyplot.xlabel('x')
    pyplot.ylabel('density')
    pyplot.title('test')
    line_ani = animation.FuncAnimation(fig1, update, fargs=(cli_args,))

    pyplot.show()

