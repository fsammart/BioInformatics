import re
from itertools import islice
import urllib.request
from exercises.aux import Aux


class Exercise4:

    @staticmethod
    def run(input_file_name, pattern, output_file_name, number_of_results):

        regex = re.compile('\|\d+\s.+\sRecName')

        output_file_name = Aux.unique_file(output_file_name)

        with open(input_file_name) as f:
            while True:
                next_n_lines = list(islice(f, number_of_results))
                chunk = ''.join(next_n_lines)
                if pattern.upper() in chunk.upper():
                    result = re.search(regex, next_n_lines[6])
                    if result:
                        full_accession = result.group(0)
                    else:
                        break
                    accession = full_accession[7:15]
                    content = urllib.request.urlretrieve("https://www.uniprot.org/uniprot/" + accession + ".fasta",
                                                         'archives/prot_sequences/'+ accession + ".fasta")
                    with open(output_file_name, "a") as out_handle:
                        print(chunk, file=out_handle)
                if not next_n_lines:
                    break
