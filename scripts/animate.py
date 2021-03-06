#!/usr/bin/env python

import argparse
import json
import numpy as np
import sys

from matplotlib import pyplot
import matplotlib.animation as animation

_color_order = [
        "steelblue", "green", "red", "cyan", "magenta", "yellow", "black"]

def parse_options():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '--frames', type=int, default=10)
    parser.add_argument(
            '--alpha', type=float, default=1.0,
            help="alpha transparency value for front histogram")
    parser.add_argument(
            '--yaxis', type=str, default="max",
            help="Set the yaxis to scale to the max observed value "
                 "(--yaxis=max), the log of the max observed value "
                 "(--yaxis=log), or a max of a specific value "
                 "(--yaxis=<FLOAT>)")

    return parser.parse_args()


def update(num, cli_args):
    global _color_order
    try:
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
                if cli_args.yaxis == "max":
                    weights = np.array(data['weights'])
                elif cli_args.yaxis == "log":
                    weights = np.array(
                            [np.log10(weight + 1.0) for weight in data['weights']])
                else:
                    max_val = float(cli_args.yaxis)
                    weights = np.array(
                            [min(max_val, weight)
                            for weight in data['weights']])

                boundaries = np.array(data['boundaries'])
                try:
                    widths = boundaries[1:] - boundaries[:-1]
                    max_weights = max(weights)
                    if ymax is None or max_weights > ymax:
                        ymax = max_weights
                    pyplot.fill_between(
                            boundaries.repeat(2)[1:-1], weights.repeat(2),
                            facecolor=memo["color"],
                            alpha=(cli_args.alpha *
                                   ((idx + 1.0) / len(memo["history"]))))
                except:
                    pass
        pyplot.ylim(0, ymax)
    except:
        pass
update.memos = {}

if __name__ == '__main__':
    cli_args = parse_options()

    fig1 = pyplot.figure()

    pyplot.xlabel('x')
    pyplot.ylabel('density')
    pyplot.title('test')
    line_ani = animation.FuncAnimation(fig1, update, fargs=(cli_args,))

    pyplot.show()

