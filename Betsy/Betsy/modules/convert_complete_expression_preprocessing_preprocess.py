from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import shutil
        in_data = antecedents
        shutil.copytree(in_data.identifier, out_path)


    def name_outfile(self, antecedents, user_options):
        return "report"



