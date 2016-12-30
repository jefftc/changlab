from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        import calc_diffexp_microarray
        
        data_node, cls_node = antecedents
        if not os.path.exists(out_path):
            os.mkdir(out_path)
        metadata = {}

        algorithm = out_attributes["de_algorithm"]
        metadata["algorithm"] = algorithm
        
        commands = calc_diffexp_microarray.calc_diffexp(
            data_node.identifier, cls_node.identifier, user_options, num_cores,
            out_path, algorithm, True)
        metadata["commands"] = commands
        metadata["num_cores"] = num_cores

        return metadata

    def name_outfile(self, antecedents, user_options):
        return "differential_expression"
