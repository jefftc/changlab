from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        from genomicode import filelib
        
        ref_dict, ref_sam, ref_bowtie1, ref_bowtie2, ref_bwa = antecedents
        
        all_ref_genomes = ref_dict, ref_sam, ref_bowtie1, ref_bowtie2, ref_bwa
        for ref in all_ref_genomes:
            # Symlink the files for each of the individual indexes.
            filelib.symlink_file_or_path_to_path(
                ref.identifier, out_path, overwrite_outpath=False)
            

    def name_outfile(self, antecedents, user_options):
        return "ref.index"
