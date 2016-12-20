from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import relabel_samples
        
        data_node, rename_node = antecedents
        metadata = {}
        cmd = relabel_samples.relabel(
            data_node.identifier, rename_node.identifier, outfile,
            user_options)
        metadata["command"] = cmd
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "control.txt"
