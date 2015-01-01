"""
Copyright: matplotlib.com example code under a PSF license (BSD compatible)
http://matplotlib.org/examples/animation/histogram.html

This example shows how to use a path patch to draw a bunch of
rectangles for an animated histogram
"""
import numpy as np
from scipy import stats

import matplotlib.pyplot as plt
plt.rcParams['animation.ffmpeg_path'] = '/home/logan/Installs/bin/ffmpeg'
import matplotlib.patches as patches
import matplotlib.path as path
import matplotlib.animation as animation
import roll
import decay_equations
import datasets

history_xlim_lower = []
history_xlim_upper = []
history_ylim_upper = []
def get_animation_function(histograms, value_generator):
    def animation_func(i):
        global history_xlim_lower, history_xlim_upper, history_ylim_upper
        plt.cla()
        # histogram our data with numpy
        colors = ["blue", "darkgreen", "magenta", "red", "yellow", "cyan"]
        for val in value_generator.next():
            for histogram in histograms:
                histogram.update(val)

        for idx, hist in enumerate(histograms):
            n, bins = hist.get_histogram()
            # In order to draw X bins, we need X + 1 bin edges. To get this,
            # we'll take the width of the first full and extend that to the
            # left of bin[0].
            bins = np.insert(bins, 0, np.array([2.0 * bins[0] - bins[1]]))

            # get the corners of the rectangles for the histogram
            left = np.array(bins[:-1])
            right = np.array(bins[1:])
            bottom = np.zeros(len(left))
            top = bottom + n
            nrects = len(left)

            # here comes the tricky part -- we have to set up the vertex and path
            # codes arrays using moveto, lineto and closepoly

            # for each rect: 1 for the MOVETO, 3 for the LINETO, 1 for the
            # CLOSEPOLY; the vert for the closepoly is ignored but we still need
            # it to keep the codes aligned with the vertices
            nverts = nrects*(1+3+1)
            verts = np.zeros((nverts, 2))
            codes = np.ones(nverts, int) * path.Path.LINETO
            codes[0::5] = path.Path.MOVETO
            codes[4::5] = path.Path.CLOSEPOLY
            verts[0::5,0] = left
            verts[0::5,1] = bottom
            verts[1::5,0] = left
            verts[1::5,1] = top
            verts[2::5,0] = right
            verts[2::5,1] = top
            verts[3::5,0] = right
            verts[3::5,1] = bottom

            barpath = path.Path(verts, codes)
            patch = patches.PathPatch(
                    barpath, facecolor=colors[idx], edgecolor='black', alpha=0.5)
            ax.add_patch(patch)

            #ax.set_xlim(-3, 3)
            #ax.set_ylim(bottom.min(), 0.03)
            history_xlim_lower.append(left[0])
            history_xlim_upper.append(right[-1])
            history_ylim_upper.append(top.max())
            xlim_lower = min(history_xlim_lower)
            xlim_upper = max(history_xlim_upper)
            ylim_upper = max(history_ylim_upper)
            if len(history_xlim_lower) > 50:
                history_xlim_lower = history_xlim_lower[1:]

            if len(history_xlim_upper) > 50:
                history_xlim_upper = history_xlim_upper[1:]

            if len(history_ylim_upper) > 50:
                history_ylim_upper = history_ylim_upper[1:]
            #ax.set_xlim(-10, 10)
            ax.set_xlim(xlim_lower, xlim_upper)
            #ax.set_ylim(bottom.min(), 0.05)
            ax.set_ylim(bottom.min(), ylim_upper)
    return animation_func


def get_histograms():
    vals = []

    alpha_count_slow = 0.001
    alpha_count_fast = decay_equations.find_optimal_decay(alpha_count_slow)
    rh_slow = roll.create_rolling_histogram_class(
            Bucket=roll.create_bucket_class(
                alpha_mu=0.01, alpha_count=alpha_count_slow),
            target_buckets=50)()
    rh_fast = roll.create_rolling_histogram_class(
            Bucket=roll.create_bucket_class(
                alpha_mu=0.01, alpha_count=alpha_count_fast),
            target_buckets=40)()
    hists = [rh_slow, rh_fast]
    return hists


def value_generator():
    dists = [stats.expon(-2, 2), stats.norm(0, 1), stats.norm(-1, 3)]
    cycles = 100
    samples_per_setting = 1
    while True:
        for _ in range(cycles):
            yield roll.gen_value(dists[0:1], samples_per_setting)
        for _ in range(cycles):
            yield roll.gen_value(dists[1:2], samples_per_setting)
        for _ in range(cycles):
            yield roll.gen_value(dists[2:3], samples_per_setting)

def foo_generator():
    for val in datasets.dataset1():
        wave = []
        for _ in xrange(10):
            wave.append(val)
        yield wave

fig, ax = plt.subplots()
animation_function = get_animation_function(
        get_histograms(), value_generator())
        #get_histograms(), foo_generator())
#FF_writer = animation.FFMpegWriter()
ani = animation.FuncAnimation(fig, animation_function, frames=None, repeat=False)
#ani.save('animated_histogram.mp4', writer=FF_writer, fps=30)
plt.show()

