from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import subprocess
        from Betsy import read_label_file
        from Betsy import module_utils
        from genomicode import config
        data_node, cls_node = antecedents
        if data_node and cls_node:
            result, label_line, second_line = read_label_file.read(
                cls_node.identifier)
            assert len(
                result) >= 2, 'for combat,there should be equal or larger than 2 classes'
            combat_path = config.combatnorm
            combat_BIN = module_utils.which(combat_path)
            assert combat_BIN, 'cannot find the %s' % combat_path
            command = ['python', combat_BIN, '-f', data_node.identifier, '-o',
                       outfile, '-label', cls_node.identifier]
            process = subprocess.Popen(command,
                                       shell=False,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            error_message = process.communicate()[1]
            if error_message:
                raise ValueError(error_message)
            assert module_utils.exists_nz(outfile), (
                'the output file %s for combat fails' % outfile
            )
        
        return False


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, cls_node = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        filename = 'signal_combat_' + original_file + '.tdf'
        return filename

    def set_out_attributes(self, antecedents, out_attributes):
        data_node, cls_node = antecedents
        new_parameters = data_node.data.attributes.copy()
        new_parameters['combat_norm'] = 'yes'
        return new_parameters


    
