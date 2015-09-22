import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as stats

# Use Helvetica as the default font family
mpl.rcParams['font.family'] = 'Helvetica2'
mpl.rcParams['font.size'] = 12

fusion_type_color = {"D":  "#AC4142",
                     "TD": "#6A9FB5",
                     "HH": "#90A959",
                     "TT": "#AA759F"}.get

"""
Test A - Clustering of breakpoints (exponential QQ plot)
Test B - Oscillation among a few CN states
Test C - Interspersed loss and retention of heterozygosity
Test D - Rearrangements affecting only one haplotype
Test E1 - Randomness of fragment joins
Test E2 - Randomness of fragment order
Test F - Ability to walk the derivative chromosome.
"""

## Test E1 - randomness of fragment joins. ##

def fusion_type_counts(fusions):
    """Takes a list of fusions, returns a pd.Series
       of counts of fusion types."""
    fusion_types = pd.Series(map(lambda x: x.type(), fusions))
    counts = fusion_types.value_counts(sort = False)
    # Put fusion types in the right order..
    ordered = counts[['D','TD','HH','TT']]
    return ordered

def plot_counts(counts, outfile, xlabel = None):
    """From a pd.Series of counts of fusion types,
       plots a histogram."""
    fig, axes = plt.subplots()
    fig.set_figwidth(3), fig.set_figheight(3)

    # Ticks face outwards
    axes.get_yaxis().set_tick_params(direction='out')

    # Only draw right and bottom ticks
    axes.yaxis.set_ticks_position('left')
    axes.xaxis.set_ticks_position('none')

    plt.locator_params(nbins = 4)

    # Only draw right and bottom spines
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['bottom'].set_visible(True)

    xticks = np.arange(4)
    width = .8

    axes.set_xlim(-(1 - width), max(xticks) + 1)

    axes.bar(left = xticks,
             width = width,
             height = counts,
             color = list(counts.index.map(fusion_type_color)),
             )

    axes.set_xticks(xticks + width/2)
    axes.set_xticklabels(counts.index,
                         size = 14,
                         ha = 'center')

    if xlabel != None:
        axes.set_xlabel(xlabel, fontsize = 15)

    plt.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)

def chisq_test(counts):
    """From a pd.Series of counts of fusion types,
       conducts a chi-squared test."""
    count_array = np.array(counts)
    chisq, p = stats.chisquare(count_array)
    return chisq, p

def test_E1(fusions, outfile, label = None):
    """Given a list of fusions, plots a histogram of counts
       of fusion types, conducts a chi-square test, and
       prints the p-value below, with an optional label below that."""
    counts = fusion_type_counts(fusions)
    chisq, p = chisq_test(counts)
    xlabel = "P = %.4f" % p
    if label != None:
        xlabel += "\n" + label
    plot_counts(counts, outfile, xlabel)


## Test F - Ability to "walk" the derivative chromosome. ##

def acc_alternating_runs(segments, walk): # Names could be improved here.
    """Recursive function that takes an empty list and a H/T
       string and returns a list of alternating segments."""
    if len(walk) == 0:
        return segments
    else:
        next_letter = walk[0]

    def conjoin(segments, next_letter):
        if segments == []:
            return [next_letter]
        else:
            last_segment = segments[-1]
            last_letter = last_segment[-1]
            if last_letter == next_letter:
                return segments + [next_letter]
            elif last_letter != next_letter:
                return segments[:-1] + [last_segment + next_letter]

    if len(walk) == 1:
        return conjoin(segments, next_letter)
    else:
        tail = walk[1:]
        return acc_alternating_runs(
                    conjoin(segments, next_letter),
                    tail)

def alternating_runs(walk):
    return acc_alternating_runs([], walk)

def modified_wald_wolfowitz(alternating_runs, N_1, N_2):
    """Returns the mean and variance of the approximate sampling
       distribution of the number of alternating runs in a sequence with
       N_1 heads and N_2 tails, as well as a one-sided p-value for
       alternating_runs. See ww_test.md for details."""

    N = N_1 + N_2
    mean_runs = 1 + (2 * N_1 * N_2 / N)
    var_runs = 2 * N_1 * N_2 * (2 * N_1 * N_2 - N) / (N**2 * (N - 1))

    # If there are R
    mean_alternating = N - mean_runs + 1
    var_alternating = var_runs

    dist = stats.norm(loc = mean_alternating,
                      scale = np.sqrt(var_alternating))
    pvalue = dist.cdf(alternating_runs)

    return mean_alternating, var_alternating, pvalue

def F_walk(walk):
    """Breaks-up a walk, prints it (and its broken-up counterpart), and
       conducts a modified Wolfowitz-Wald test."""
    head_prop = (walk.count("H") / float(len(walk)))

    alt_runs = alternating_runs(walk)
    N_1 = walk.count("H")
    N_2 = walk.count("T")
    N = N_1 + N_2
    average_run_length = N / float(len(alt_runs))

    mean, var, pvalue = modified_wald_wolfowitz(len(alt_runs), N_1, N_2)

    print walk
    print "|".join([r for r in alt_runs])
    print "heads: %s; tails: %s" % (N_1, N_2)
    print ("alternating runs: %s; average run length: %.2f"
            % (len(alt_runs), average_run_length))
    print "expected alternating runs: ~ %s; sd: %s" % (mean, np.sqrt(var))
    print "p-value: %.4g\n" % pvalue
