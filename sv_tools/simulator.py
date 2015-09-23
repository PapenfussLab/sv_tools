import numpy as np

import sv_data
import sv_diagram as sv_d

def map_kmers(f, k):
    """ Takes a list function f and returns a function that applies
        f to k-mers of a list, returning the results as a list with
        None values discarded.
    """
    def g(input_list, *args, **kwargs):
        outputs = [f(input_list[i:i+k], *args, **kwargs)
                   for i in range(len(input_list) + 1 - k)]
        return [x for x in outputs if x != None]

    return g

def map_kbins(f, k):
    def g(input_list, *args, **kwargs):
        # Length of the input list must be a multiple of k.
        assert len(input_list) % k == 0
        outputs = [f(input_list[i:i+k], *args, **kwargs)
                   for i in range(0, len(input_list), k)]
        return [x for x in outputs if x != None]

    return g

# Sequence of letters to rearranged positions

def pair_to_letters(pair):
    letter, possible_tick = pair

    if not letter.isalpha():
        return None
    elif possible_tick != "'":
        return letter
    else:
        return letter + possible_tick

def letters_to_letterlist(letters):
    letters += "A"
    to_letterlist = map_kmers(pair_to_letters, 2)
    return to_letterlist(letters)

def pair_to_positions(pair, length = 10):
    letter, possible_tick = pair
    positions = list(np.arange(length) + (ord(letter) * length))
    inverted = list(positions[::-1])

    if not letter.isalpha():
        return None
    elif possible_tick != "'":
        return positions
    else:
        return inverted

def letters_to_positions(letters):
    letters += "A" # Not read.
    pairs_to_positions = map_kmers(pair_to_positions, 2)
    positions = [x for letter_sequence in pairs_to_positions(letters)
                   for x in letter_sequence]
    return positions

# Rearranged positions to sequence of letters

def positions_to_letter(positions, length = 10):
    assert len(positions) == length
    base_position = min(positions[0], positions[-1])
    if base_position == positions[0]:
        inverted = False
    elif base_position == positions[-1]:
        inverted = True
    else:
        print positions[0], positions[-1], base_position

    letter = chr(base_position / len(positions))

    if not inverted:
        return letter
    else:
        return letter + "'"

def positions_to_letters(positions):
    positions_to_list = map_kbins(positions_to_letter, 10)
    list_of_letters = positions_to_list(positions)
    return "".join(list_of_letters)

positions_to_ticks = map_kbins(np.mean, 10)

# Fusions from rearranged chromosome.

def detect_fusions(sites):
    """ Takes a list of four sites, and returns either None, or a
        fusion-tuple based on a length-four paired-end read, e.g.

        01[2398]7
           -><-
           T  T
    """
    assert len(sites) == 4

    breakdiff = abs(sites[1] - sites[2])
    diff1 = sites[0] - sites[1]
    diff2 = sites[2] - sites[3]

    if breakdiff == 1:
        return None
    else:
        # Differences should be 1 or -1 normally.
        strand1 = {-1:"+", 1:"-"}.get(diff1, "?")
        strand2 = {1:"+", -1:"-"}.get(diff2, "?")
        bp1 = sv_data.Breakpoint(chrom = "",
                              pos = sites[1] * 1e6,
                              strand = strand1)
        bp2 = sv_data.Breakpoint(chrom = "",
                              pos = sites[2] * 1e6,
                              strand = strand2)
        return sv_data.Fusion(bp1, bp2)

get_fusions = map_kmers(detect_fusions, 4)

# Copy number from rearranged chromosome

def get_x_cn(positions):
    counts = [(p, positions.count(p)) for p in positions]
    x_tuple, cn_tuple = zip(*counts)
    x, cn = list(x_tuple), list(cn_tuple)
    return x, cn

## Campbellgrams ##

def simulate_sv_diagram(
        letters, outfile = None,
        **kwargs):

    if outfile == None:
        outfile = "../output/simulation/simulation_%s.pdf" % letters

    ### Simulation-specific stuff
    positions = letters_to_positions(letters)
    fusions = get_fusions(positions)
    x, cn = get_x_cn(positions)
    kwargs['yticks'] = range(max(cn) + 2)
    kwargs['ymax'] = max(cn) + 1
    kwargs['ymin'] = 0
    kwargs['xlabel'] = letters
    ###

    fig = sv_d.setup_figure()
    cn_axes, fusion_axes = sv_d.sv_diagram_axes()

    # Copy number

    sv_d.plot_cn(cn_axes, x, cn)

    sv_d.set_cn_axes_options(cn_axes, x, cn, kwargs)
    sv_d.set_cn_axes_aesthetics(cn_axes)
    sv_d.plt.minorticks_off()


    ### Simulation-specific stuff
    x_range = range(min(x), max(x) + 1)
    x_letters = letters_to_letterlist(positions_to_letters(x_range))
    x_ticks = positions_to_ticks(x_range)
    cn_axes.set_xticks(x_ticks, minor = True)
    cn_axes.set_xticklabels(x_letters, minor = True)

    sv_d.plt.setp(cn_axes.get_xticklabels(), visible=False)

    ###

    # Fusions

    sv_d.setup_fusion_axes(fusion_axes, min(x), max(x))
    for fusion in fusions:
        sv_d.plot_fusion(cn_axes, fusion_axes, fusion)

    # Ensure everything fits
    sv_d.plt.tight_layout()

    # Output

    fig.savefig(outfile)
    sv_d.plt.close(fig)
