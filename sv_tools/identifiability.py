import simulator as sim
import itertools as it
import numpy     as np
import networkx  as nx
import re

# Helper functions

def parse_string(string):
    return re.findall(".\'?", string)

def flip(token):
    if len(token) == 2:
        return token[0]
    elif len(token) == 1:
        return token + "'"

def reverse_string(string):
    reversed_tokens = map(flip, reversed(parse_string(string)))
    return "".join(reversed_tokens)

# Main class

class ChromString(object):
    """ This "equivalence class" allows the use of sets to add speed 
        and eliminate redundancy. Currently only using it within 
        functions.
        
        Create using a string like "AB'CE". """
    
    def __init__(self, raw_string):
        if raw_string.count("'") < \
           reverse_string(raw_string).count("'"):
            self.string = raw_string
        elif raw_string.count("'") > \
             reverse_string(raw_string).count("'"):
            self.string = reverse_string(raw_string)
        elif raw_string <= reverse_string(raw_string):
            self.string = raw_string
        else:
            self.string = reverse_string(raw_string)

    @staticmethod
    def from_tokens(tokens):
        return ChromString("".join(tokens))
    
    # Main methods -- all non-mutating I hope.
    
    def tokens(self):
        return parse_string(self.string)
    
    def reversed_string(self):
        return reverse_string(self.string)
        
    def flipped_indices(self, indices):
        tokens = self.tokens()
    
        def flip_if(index, token):
            if index in indices:
                return flip(token)
            else:
                return token
    
        new_tokens = [flip_if(i, t) for i, t in enumerate(tokens)]
        
        return ChromString("".join(new_tokens))
        
    def deleted_indices(self, indices):
        tokens = self.tokens()
        
        new_tokens = [t for i, t in enumerate(tokens)
                      if not (i in indices)]
                      
        return ChromString.from_tokens(new_tokens)
        
    def plot_sv_diagram(self, outfile = None):
        sim.simulate_sv_diagram(self.string, outfile = None)
    
    def __repr__(self):
        return self.string
        
    def __eq__(self, other):
        return self.string == other.string
        
    def __hash__(self):
        return hash(self.__repr__())

# Helper function(s)

def powerset(iterable):
    """ powerset([1,2,3]) --> 
        () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3) """
    s = list(iterable)
    return it.chain.from_iterable(it.combinations(s, r) 
                                  for r in range(len(s) + 1))

# --- Permutations, inversions, deletions, and rearrangements

def all_permutations(chrom_string):   
    "Returns all permutations of the characters of a string."
    permuted_tokens = it.permutations(chrom_string.tokens())
    return set(map(ChromString.from_tokens, permuted_tokens))

def all_inversions(chrom_strings):
    out_strings = (s.flipped_indices(indices)
                   for s in chrom_strings
                   for indices in powerset(range(len(s.tokens()))))
    
    return set(out_strings)

def all_deletions(chrom_strings):
    out_strings = (s.deleted_indices(indices)
                   for s in chrom_strings
                   for indices in powerset(range(len(s.tokens()))))
    
    return set(out_strings)
    
def all_rearrangements(chrom_string):
    # This is highly inefficient: need to avoid redundancy.
    return set(
           all_inversions(all_deletions(all_permutations(chrom_string)))
           )
    
# --- Identifiability

def fusion_set(chrom_string):
    return set(sim.get_fusions(sim.letters_to_positions(chrom_string.string)))
    
def sv_diagram_data(chrom_string):
    if len(chrom_string.string) > 0:
        x, cn = sim.get_x_cn(sim.letters_to_positions(chrom_string.string))
        letters = set([l for l in chrom_string.string if l.isalnum()])
        return letters, cn, fusion_set(chrom_string)
    else:
        return None

def find_clashes(chrom_strings, mapping_function):
    d = {s:mapping_function(s) for s in chrom_strings}
    clash_pairs = ((x,y) for x in d for y in d 
                         if (x != y and d[x] == d[y]))
    clash_components = nx.connected_components(nx.Graph(clash_pairs))
    return clash_components
    