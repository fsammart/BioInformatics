import sys
import os
import itertools as itertools

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

import regex as re

INPUT_PARAMS_COUNT = 3

#returns the versioned name in case filename already existed.
def unique_file(name):
	basename, ext = os.path.splitext(name)
	actualname = "%s%s" % (basename, ext)
	c = itertools.count()
	next(c)
	while os.path.exists(actualname):
		actualname = "%s (%d)%s" % (basename, next(c), ext)

	return actualname

#Processing starts
if len(sys.argv) != INPUT_PARAMS_COUNT:
	print("Invalid params: 1) Input file path 2) Output file path")
	sys.exit(1)

input_file_name = sys.argv[1]
# if output file name already exists we append a number '(n)' indicating version
output_file_name = unique_file(sys.argv[2])

# we are analyzing just one gene for this disease, but the algorithm can process multiple genes.
records = SeqIO.parse(input_file_name, "gb")

for record in records:
	record_str = str(record.seq)
	rev_record_str = str(record.seq.reverse_complement())
	startCodon = re.compile('ATG')
	nuc_chain = record_str.replace('\n','')
	rev_nuc_chain = rev_record_str.replace('\n','')
	longest = (0,)
	for seq in nuc_chain, rev_nuc_chain:
		for m in startCodon.finditer(seq, overlapped=True):
			rec = Seq(seq)
			if len(rec[m.start():].translate(to_stop=True)) > longest[0]:
				#secuence found is longer than the stored one
				pro = rec[m.start():].translate(to_stop=True)
				longest = (len(pro),str(pro))

	protein_record = SeqRecord(Seq(longest[1], IUPAC.protein), id=record.id, description= record.description + " protein translation")

	with open(output_file_name, "a") as output_handle:
		SeqIO.write(protein_record, output_handle, "fasta")