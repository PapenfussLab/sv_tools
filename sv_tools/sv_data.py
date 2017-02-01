import pandas as pd

### Main classes ###

class Breakpoint(object):
    """A Breakpoint is essentially a chromosomal coordinate with an
       orientation."""
    def __init__(self, chrom, pos, strand):
        self.chrom = chrom
        self.pos = pos
        self.strand = strand

    def orientation(self):
        """Orientation, as determined by strand."""
        return {'+': 'T', '-': 'H'}.get(self.strand)

    def pos_scaled(self):
        """Position in Mb."""
        return self.pos / 1e6

    def __repr__(self):
        return "%s:%s(%s)" % (self.chrom,
                              str(self.pos),
                              self.orientation())

    # Allows "=="
    def __eq__(self, other):
        return self.__dict__ == other.__dict__
    # Allows set() to determine equality.
    def __hash__(self):
            return hash(self.__repr__())

class Fusion(object):
    """A Fusion is a pair of breakpoints."""
    def __init__(self, first_bp, second_bp):
        ordered = sorted([first_bp, second_bp],
                         key = lambda x: x.pos)
        self.bp1 = ordered[0]
        self.bp2 = ordered[1]
        # self.reads = reads  # Not used
        # self.gap = gap      # Not used

    def orientations(self):
        """Fusion type in TH/HT/HH/TT form."""
        return self.bp1.orientation() + self.bp2.orientation()

    def type(self):
        """Fusion type in D/TD/HH/TT form."""
        return {"TH": "D",
                "HT": "TD",
                "HH": "HH",
                "TT": "TT"}.get(self.orientations())

    def __repr__(self):
        # Ugly!
        return self.bp1.__repr__() + "--->" + self.bp2.__repr__()

    # Allows "=="
    def __eq__(self, other):
        return self.__dict__ == other.__dict__
    # Allows set() to determine equality.
    def __hash__(self):
            return hash(self.__repr__())


def df_from_txt(txt_file):
    """Returns a Pandas data frame from tab-delimited text"""

    return pd.read_csv(txt_file,
                       sep = "\t",
                       header = 0)

### Getting Fusions and Breakpoints from files. Key interfaces. ###

def fusion_from_row(row):
    """Defines the interface between rows of
       the data[frame], and Breakpoint/Fusion objects."""

    bp1 = Breakpoint(row['chrom1'], row['pos1'], row['strand1'])
    bp2 = Breakpoint(row['chrom2'], row['pos2'], row['strand2'])
    return Fusion(bp1, bp2)

def fusions(fusion_data):
    """Returns a list of Fusions from a Pandas dataframe of
       fusions, possibly restricting attention to one chromosome."""

    fusion_series = fusion_data.apply(fusion_from_row, axis = 1)
    return list(fusion_series)

def breakpoints(fusion_data):
    """Returns a list of Breakpoints from a Pandas dataframe
       of fusions."""

    return [bp for fusion in fusions(fusion_data)
               for bp in (fusion.bp1, fusion.bp2)]

# Convenience functions

def get_fusions(filename, chrom):
    """Get the fusions for a given chromosome from a file."""
    fs = fusions(df_from_txt(filename))
    return [f for f in fs if f.bp1.chrom == chrom and
                             f.bp2.chrom == chrom]

def get_breakpoints(filename, chrom):
    """Get the breakpoints for a given chromosome from a file."""
    all_breaks = breakpoints(df_from_txt(filename))
    breaks = [bp for bp in all_breaks if bp.chrom == chrom]
    sorted_breaks = sorted(breaks, key = lambda x: x.pos)
    return sorted_breaks

### Getting copy number data

def get_x_cn(filename, chrom):
    """Get the copy number info for a given chromosome from a file."""
    # Rewritten for .bed files
    fields = ["chrom", "start", "end", "name", "CN"]
    df = pd.read_csv(filename,
                     sep = '\t',
                     names = fields)
    df_chrom = df[df['chrom'] == chrom]
    x = df_chrom['start']
    depth = df_chrom['CN']
    return x, depth
