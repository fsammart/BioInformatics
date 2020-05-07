from Bio.Align.Applications import ClustalwCommandline
from exercises.aux import Aux


class Exercise3:

    @staticmethod
    def run(input_file_name, output_file_name):
        output_file_name = Aux.unique_file(output_file_name)
        print('Saving result in ' + output_file_name + ' file')
        clustalw_cline = ClustalwCommandline("clustalo", infile=input_file_name, outfile=output_file_name)
        stdout, stderr = clustalw_cline()
