import os
from ast import literal_eval

path_to_breakpoints = '/Users/franciscosammartino/ITBA/BioInformatics/Sorting-by-Fragmentation-Weighted-Rearrangements/'


class BreakpointsReversal:
    @staticmethod
    def run(disordered_array):
        n = len(literal_eval(disordered_array))
        disordered_array = disordered_array.replace('[','').replace(']','')

        os.system(path_to_breakpoints + './prog_o -a sr_g -v4 -n' + str(n) + ' -p '+ disordered_array)
