import os

def reversal(seq):
    seq.reverse()
    seq = map(lambda x: -x, seq)
    return seq


def greedy_reordering(lst):
    reversals = []
    steps = []

    for i in range(len(lst)):
        # We check if the element is in the right position.
        if lst[i] != i+1:
            if i+1 in lst:
                # We search for the index of the element (i+1).
                id_i = lst.index(abs(i+1))
            else:
                # If we cannot find it, means it has the negative sign.
                id_i = lst.index(-abs(i+1))
            # We reverse from the current position, to the one with the element.
            steps.append([i,id_i])
            lst[i: id_i + 1] = reversal(lst[i: id_i + 1])
            reversals.append(lst[:])
        # We check if the element has the right sign.
        if lst[i] == -(i+1):
            lst[i] = i+1
            reversals.append(lst[:])
            steps.append([i,i])

    return reversals,steps

def count_breakpoints(lst):
    '''
    Returns number of breakpoint in sequence.
    (pi,pi+1) is a breakpoint <-> pi+1 − pi is equal to 1. Ej.: (4,5) or
    (-5,-4)
    '''
    count = 0;

    for i in range(len(lst) - 1):
        if (lst[i+1] - lst[i] != 1):
            count+= 1

    return count


path_to_grappa = '/Users/franciscosammartino/ITBA/BioInformatics/GRAPPA20/'
class Greedy:
    @staticmethod
    def run(disordered_array):
        s = '>S1\n'
        base = '\n>S2\n'
        for i,n in enumerate(disordered_array):
            s += str(n) + ' '
            base += str(i+1) + ' '
        print('Breakpoints = ' + str(count_breakpoints(disordered_array)))
        reversals,steps = greedy_reordering(disordered_array)
        print('Reversals Count = ' + str(len(reversals)))
        print('\n')
        for idx,elem in enumerate(reversals):
            print('Reversal Idx:',steps[idx],'-> ', elem)
        os.system('echo "' + s + base + '" >' + path_to_grappa + 'modified.so' )
        print('\n\nMinimum Reversals Bader, D.A., Moret, B.M.E., and Yan, M. \n')
        os.system(path_to_grappa + './invdist -f ' + path_to_grappa + 'modified.so | grep "Inversion Distance"')
        print('\n')
