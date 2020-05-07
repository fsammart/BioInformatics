from Bio.Align.Applications import ClustalOmegaCommandline
from exercises.aux import Aux


class Exercise3:

    @staticmethod
    def run(input_file_name, output_file_name):

        output_file_name = Aux.unique_file(output_file_name)
        print(input_file_name)
        print(output_file_name)
        clustalw_cline = ClustalOmegaCommandline(infile=input_file_name, outfile=output_file_name)
        stdout, stderr = clustalw_cline()