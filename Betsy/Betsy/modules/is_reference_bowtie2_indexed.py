from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        from genomicode import alignlib

        alignlib.standardize_reference_genome(
            in_data.identifier, out_path, use_symlinks=True)


    def name_outfile(self, antecedents, user_options):
        return "reference.bowtie2"

    
    def set_out_attributes(self, in_data, out_attributes):
        from genomicode import alignlib

        ref = alignlib.create_reference_genome(in_data.identifier)
        is_indexed = "yes"
        if not ref.bowtie2_indexes:
            is_indexed = "no"
        
        attrs = out_attributes.copy()
        attrs["bowtie2_indexed"] = is_indexed
        return attrs
        

