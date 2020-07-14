import os


path_to_circular = '/Users/franciscosammartino/ITBA/BioInformatics/Sorting-by-Fragmentation-Weighted-Rearrangements/adapted-algorithms/GRIMM-2.01/'
class CircularReversal:
    @staticmethod
    def run(disordered_array):
        s = '>CHR1\n'
        base = '\n>CHR2\n'
        for i,n in enumerate(disordered_array):
            s += str(n) + ' '
            base += str(i+1) + ' '
        s += '$'
        os.system('echo "' + s + base + '" >' + path_to_circular + 'file.data' )
        os.system(path_to_circular + './grimm -C -f ' + path_to_circular + 'file.data')
        print('\n')
