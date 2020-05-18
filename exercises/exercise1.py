from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from exercises.aux import Aux

import regex as re

CODON_LENGTH = 3


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
                    rec = Seq(seq)[m.start():]
                    # we check if sequence is multiple of 3 (codon length)
                    if len(rec) % CODON_LENGTH == 0:
                        aux = rec.translate(to_stop=True)
                        if len(aux) > longest[0]:
                            pro = aux
                            longest = (len(pro), str(pro))

            protein_record = SeqRecord(Seq(longest[1], IUPAC.protein), id=record.id,
                                       description=record.description + " protein translation")

            Aux.save_file(output_file_name, protein_record)
