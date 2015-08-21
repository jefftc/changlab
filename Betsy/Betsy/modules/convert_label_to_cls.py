from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import arrayio
        from Betsy import read_label_file
        from Betsy import module_utils
        data_node, cls_node = antecedents
        if data_node and cls_node:
            f = file(cls_node.identifier, 'rU')
            text = f.readlines()
            f.close()
            text = [i.rstrip() for i in text]
            label_dict = {}
            for line in text:
                words = line.split('\t')
                if words[1] in label_dict:
                    label_dict[words[1]].append(words[0])
                else:
                    label_dict[words[1]] = [words[0]]
            class_names = label_dict.keys()
            M = arrayio.read(data_node.identifier)
            column_names = M.col_names('_SAMPLE_NAME')
            label_line = [0] * len(column_names)
            for i in range(len(class_names)):
                sample_names = label_dict[class_names[i]]
                for sample_name in sample_names:
                    index = column_names.index(sample_name)
                    label_line[index] = str(i)
            read_label_file.write(outfile, class_names, label_line)
            assert module_utils.exists_nz(outfile), (
                'the output file %s for convert_label_to_cls fails' % outfile
            )
        
        return False


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, cls_node = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        filename = 'class_label_' + original_file + '.cls'
        return filename


    def hash_input(self, pipeline, antecedents, out_attributes, user_options):
        from Betsy import module_utils
        data_node, cls_node = antecedents
        identifier = data_node.identifier
        return module_utils.hash_input(identifier, pipeline, out_attributes,
                                             user_options)
