import errno
import itertools
import os

from Bio import SeqIO


class Aux:

    @staticmethod
    # returns the versioned name in case filename already existed.
    def unique_file(name):
        if not os.path.exists(os.path.dirname(name)):
            try:
                os.makedirs(os.path.dirname(name))
            except OSError as exc:  # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise
        basename, ext = os.path.splitext(name)
        current_name = "%s%s" % (basename, ext)
        c = itertools.count()
        next(c)
        while os.path.exists(current_name):
            current_name = "%s_%d%s" % (basename, next(c), ext)

        return current_name

    @staticmethod
    def save_file(output_file_name, protein_record):
        # if output file name already exists we append a number '(n)' indicating version
        output_file_name = Aux.unique_file(output_file_name)

        with open(output_file_name, "a") as output_handle:
            SeqIO.write(protein_record, output_handle, "fasta")
