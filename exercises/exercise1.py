import sys
import os
import itertools as itertools
import errno

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

import regex as re


# returns the versioned name in case filename already existed.
def unique_file(name):
    basename, ext = os.path.splitext(name)
    actualname = "%s%s" % (basename, ext)
    c = itertools.count()
    next(c)
    while os.path.exists(actualname):
        actualname = "%s (%d)%s" % (basename, next(c), ext)

    return actualname

def save_file(output_file_name, protein_record):
    # if output file name already exists we append a number '(n)' indicating version
    output_file_name = unique_file(output_file_name)
    
    if not os.path.exists(os.path.dirname(output_file_name)):
        try:
            os.makedirs(os.path.dirname(output_file_name))
        except OSError as exc:  # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    with open(output_file_name, "a") as output_handle:
        SeqIO.write(protein_record, output_handle, "fasta")


class Exercise1:

    @staticmethod
    def run(input_file_name, output_file_name):

        # we are analyzing just one gene for this disease, but the algorithm can process multiple genes.
        records = SeqIO.parse(input_file_name, "gb")

        for record in records:
            record_str = str(record.seq)
            rev_record_str = str(record.seq.reverse_complement())
            startCodon = re.compile('ATG')
            nuc_chain = record_str.replace('\n', '')
            rev_nuc_chain = rev_record_str.replace('\n', '')
            longest = (0,)
            for seq in nuc_chain, rev_nuc_chain:
                for m in startCodon.finditer(seq, overlapped=True):
                    rec = Seq(seq)
                    if len(rec[m.start():].translate(to_stop=True)) > longest[0]:
                        # secuence found is longer than the stored one
                        pro = rec[m.start():].translate(to_stop=True)
                        longest = (len(pro), str(pro))

            protein_record = SeqRecord(Seq(longest[1], IUPAC.protein), id=record.id,
                                       description=record.description + " protein translation")

            save_file(output_file_name, protein_record)
