import sys
import os
from Bio import SeqIO


def get_motifs(fasta_file, output_file):
    with(open(fasta_file)) as file_content:
        f = open('temp_fasta.fasta', 'w')
        f.write(file_content.read())
        f.close()

        os.system('prosextract -prositedir archives/prosite')
        os.system('patmatmotifs temp_fasta.fasta ' + output_file)
        os.remove('temp_fasta.fasta')


def get_orf(gbk, output_file_orf):
    SeqIO.convert(gbk, "genbank", "temp.fasta", "fasta")
    os.system("getorf temp.fasta " + output_file_orf)
    os.remove("temp.fasta")


class Exercise5:

    @staticmethod
    def run(fasta_file, gb_file, output_file):
        if fasta_file is not None:
            get_motifs(fasta_file, output_file)
        else:
            get_orf(gb_file, output_file)
