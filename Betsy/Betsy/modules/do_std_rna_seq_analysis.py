from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        bam_node, \
                  fastqc_summary1_node, fastqc_folder1_node, \
                  fastqc_summary2_node, fastqc_folder2_node, \
                  rseqc_node, \
                  signal1_node, \
                  aligned_reads_node, \
                  signal2_node, \
                  htseq_reads_node = antecedents
        filelib.safe_mkdir(out_path)

        FILES = [
            (bam_node.identifier, False, "alignment.bam"),
            (fastqc_summary1_node.identifier, True, "fastqc.no_trim.xls"),
            (fastqc_folder1_node.identifier, False, "fastqc.no_trim"),
            (fastqc_summary2_node.identifier, True, "fastqc.trim.xls"),
            (fastqc_folder2_node.identifier, False, "fastqc.trim"),
            (rseqc_node.identifier, False, "RSeQC"),
            (signal1_node.identifier, True, "expression.tpm"),
            (aligned_reads_node.identifier, True, "aligned.xls"),
            (signal2_node.identifier, True, "expression.counts"),
            (htseq_reads_node.identifier, True, "mapped.htseq.txt"),
            ]

        for x in FILES:
            orig_filename, is_file, new_file = x
            new_filename = os.path.join(out_path, new_file)

            # Copy or link the data into the right place.
            if is_file:
                filelib.assert_exists_nz(orig_filename)
            else:
                assert filelib.dir_exists(orig_filename), \
                       "Directory not found or not directory: %s" % \
                       orig_filename
            os.symlink(orig_filename, new_filename)
        
    
    def name_outfile(self, antecedents, user_options):
        return "rna_seq"


    
