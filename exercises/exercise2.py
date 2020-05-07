import os

from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
from exercises.aux import Aux

E_VALUE_THRESH = 10
RESULTS_XML = "results.xml"
PROT_DB = "archives/swissprot/swissprot"
CMD_BLAST = "/usr/local/ncbi/blast/bin/blastp"


class Exercise2:

    @staticmethod
    def run(input_file_name, output_file_name, online=False):
        fasta_string = open(input_file_name).read()

        if online:
            result_handle = NCBIWWW.qblast("blastp", "swissprot", fasta_string, expect=E_VALUE_THRESH)

            with open(RESULTS_XML, "w") as out_handle:
                out_handle.write(result_handle.read())
            result_handle.close()
            result_handle = open(RESULTS_XML)
        else:
            blastx_cline = NcbiblastxCommandline(cmd=CMD_BLAST, query=input_file_name, db=PROT_DB,
                                                 evalue=E_VALUE_THRESH, out=RESULTS_XML, outfmt=5)
            stdout, stderr = blastx_cline()
            result_handle = open(RESULTS_XML)

        blast_records = NCBIXML.parse(result_handle)

        output_file_name = Aux.unique_file(output_file_name)

        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    with open(output_file_name, "a") as f:
                        print("****Blast Result****", file=f)
                        print("sequence:", alignment.title, file=f)
                        print("length:", alignment.length, file=f)
                        print("e value:", hsp.expect, file=f)
                        print("gaps:", hsp.gaps, file=f)
                        print("identities:", hsp.identities, file=f)
                        print("positives:", hsp.positives, file=f)
                        print("score:", hsp.score, file=f)
                        print(hsp.query[0:75] + "...", file=f)
                        print(hsp.match[0:75] + "...", file=f)
                        print(hsp.sbjct[0:75] + "...", file=f)
        if os.path.exists("results.xml"):
            os.remove("results.xml")
