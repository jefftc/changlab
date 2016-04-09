from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_filename):
        from genomicode import alignlib
        from genomicode import genomelib

        ref = alignlib.create_reference_genome(in_data.identifier)
        metadata = {}

        # Format:
        # <chr>:<start>-<end>
        # coordinates are 1-based, inclusive.
        handle = open(out_filename, 'w')
        for x in genomelib.read_fasta_many(ref.fasta_file_full):
            title, sequence = x
            # chrom is just the first part of the title.
            chrom = title.split()[0]
            start = 1
            end = len(sequence)
            print >>handle, "%s:%d-%d" % (chrom, start, end)
        handle.close()
        return metadata
    
    def name_outfile(self, antecedents, user_options):
        return "whole_genome.intervals"
