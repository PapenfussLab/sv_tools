import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from matplotlib.ticker import ScalarFormatter

# Use Helvetica as the default font family
mpl.rcParams['font.family'] = 'Helvetica2'
mpl.rcParams['font.size'] = 10.5

from mpl_toolkits.axes_grid1 import ImageGrid, AxesGrid

fusion_type_color = {"D":  "#AC4142",
                     "TD": "#6A9FB5",
                     "HH": "#90A959",
                     "TT": "#AA759F"}.get


### Campbell-grams ###

def sv_diagram_axes():
    """Set up a stacked pair of axes"""
    gs = gridspec.GridSpec(2, 1,
                           height_ratios=[.8,1],
                           hspace = 0)

    gs.update(hspace = 0) # No gap between axes.
    cn_axes = plt.subplot(gs[1])
    fusion_axes = plt.subplot(gs[0], sharex = cn_axes)
    return cn_axes, fusion_axes

## CN axes ##

def set_cn_axes_aesthetics(cn_axes):
    """Set main aesthetics of the CN axes"""
    # Ticks face outwards
    cn_axes.get_yaxis().set_tick_params(direction='out', which = 'both')
    cn_axes.get_xaxis().set_tick_params(direction='out', which = 'both')

    # Only draw right and bottom spines
    cn_axes.spines['top'].set_visible(False)
    cn_axes.spines['right'].set_visible(False)

    # Only draw right and bottom ticks
    cn_axes.yaxis.set_ticks_position('left')
    cn_axes.xaxis.set_ticks_position('bottom')

    """
    # On the y axis, only draw spines between ticks. Doesn't work.
    yticks = cn_axes.get_yaxis().get_majorticklocs()
    print yticks
    cn_axes.spines['left'].set_bounds(min(yticks), max(yticks))
    """

    # No decimal places for CN.
    cn_axes.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))

    # Avoid minor ticks -- particularly for the log scale.
    cn_axes.minorticks_off()

def set_cn_axes_options(cn_axes, x, cn, kwargs):
    """Sets various options for the CN axes, including labels,
       whether to log-transform the y-axis, etc. Could add zebra
       stripe options ..."""
    # xmin, xmax, ymin, ymax, xticks, yticks, xlabel, ylabel

    if kwargs.get("logbase", None) != None:
        cn_axes.set_yscale("log", basey = kwargs["logbase"])
        cn_axes.yaxis.set_major_formatter(ScalarFormatter())
        plt.minorticks_off()

    if (kwargs.get("xmin", None) != None) and (kwargs.get("xmax", None) != None):
        cn_axes.set_xlim(kwargs["xmin"], kwargs["xmax"])
        cn_axes.spines['bottom'].set_bounds(kwargs["xmin"], kwargs["xmax"])
    else:
        cn_axes.set_xlim(min(x), max(x))
        cn_axes.spines['bottom'].set_bounds(min(x), max(x))

    if kwargs.get("ymin", None) != None:
        cn_axes.set_ylim(bottom = kwargs["ymin"])

    if kwargs.get("ymax", None) != None:
        cn_axes.set_ylim(top = kwargs["ymax"])

    if kwargs.get("xticks", None) != None:
        cn_axes.set_xticks(kwargs["xticks"])
        cn_axes.set_xticklabels(kwargs["xticks"])

    if kwargs.get("yticks", None) != None:
        yticks = kwargs["yticks"]
        cn_axes.set_yticks(yticks)
        cn_axes.set_yticklabels(yticks)
        cn_axes.spines['left'].set_bounds(min(yticks), max(yticks))

    if kwargs.get("xlabel", None) != None:
        cn_axes.set_xlabel(kwargs["xlabel"])
    else:
        cn_axes.set_xlabel("Position on chromosome ?? (Mb)")

    if kwargs.get("ylabel", None) != None:
        cn_axes.set_ylabel(kwargs["ylabel"])
    else:
        cn_axes.set_ylabel("Estimated copy number")

def plot_cn(cn_axes, x, cn):
    """Plot copy number"""
    cn_axes.plot(x, cn, 'o', markersize=1, color = 'black', alpha = .4)

## Fusion axes ##

def setup_fusion_axes(fusion_axes, xmin, xmax):
    """Sets up the fusion axes."""

    fusion_axes.set_frame_on(False)
    fusion_axes.xaxis.set_visible(False)
    fusion_axes.yaxis.set_visible(False)

    fusion_axes.set_ylim([0, 5])

    # Lines for D/TD; HH/TT.
    fusion_axes.hlines([2, 4], 0, xmax,
                       linestyles = "dashed", alpha = .3)

    offset = 2.5 # Need to fix this to use absolute coordinates
    fusion_axes.text(xmin - offset , 4 + .2, "D", ha = 'right')
    fusion_axes.text(xmin - offset , 4 - .2, "TD", ha = 'right',
                                                   va = 'top')
    fusion_axes.text(xmin - offset , 2 + .2, "HH", ha = 'right')
    fusion_axes.text(xmin - offset , 2 - .2, "TT", ha = 'right',
                                                   va = 'top')

def plot_fusion(cn_axes, fusion_axes, fusion):
    """Plots a fusion as an arc on the fusion_axis with vertical lines
       on both axes"""

    x_coords = [fusion.bp1.pos_scaled(), fusion.bp2.pos_scaled()]

    fusion_type = fusion.type()

    height = {"D": 4,
              "TD": 4,
              "HH": 2,
              "TT": 2}

    angles = {"D": (0, 180),
              "TD": (180, 0),
              "HH": (0, 180),
              "TT": (180, 0)}

    ymin, ymax = cn_axes.get_ylim()

    def plot_vlines(axes, y_lower, y_upper):
        axes.vlines(x_coords, y_lower, y_upper,
                   colors=fusion_type_color(fusion_type),
                   alpha = 0.25) # Questionable stylistic choice ...

    plot_vlines(cn_axes, ymin, ymax)
    plot_vlines(fusion_axes, 0, height[fusion_type])


    e = patches.Arc(
                xy = [np.mean(x_coords),
                      height[fusion_type]],
                width = abs(x_coords[1] - x_coords[0]),
                height = 1.5,
                theta1 = angles[fusion_type][0],
                theta2 = angles[fusion_type][1],
                color = fusion_type_color(fusion_type),
                alpha = .7)

    fusion_axes.add_artist(e)

## Putting it all together ##

def setup_figure(width = 3.5, height = 3.3, dpi = 1000):
    """Sets up the 'canvas'."""
    fig = plt.figure()
    fig.set_dpi(dpi)
    fig.set_figwidth(width)
    fig.set_figheight(height)
    return fig

def plot_sv_diagram(x, cn, fusions, outfile, **kwargs):
    """Plots a Campbell-gram with default-y settings.
       Key word arguments are aesthetic options which can be
       safely left blank:
       xmin, xmax, ymin, ymax, xticks, yticks,
       logbase, xlabel, ylabel
    """

    # X axis is in Mb
    x = x / 1e6

    fig = setup_figure()
    cn_axes, fusion_axes = sv_diagram_axes()

    # Copy number

    plot_cn(cn_axes, x, cn)

    set_cn_axes_options(cn_axes, x, cn, kwargs)
    set_cn_axes_aesthetics(cn_axes)
    plt.minorticks_off()

    # Fusions

    setup_fusion_axes(fusion_axes, min(x), max(x))
    for fusion in fusions:
        plot_fusion(cn_axes, fusion_axes, fusion)

    # Ensure everything fits
    plt.tight_layout()

    # Output

    fig.savefig(outfile)
    plt.close(fig)
