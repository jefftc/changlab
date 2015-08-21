from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import math
        from Betsy import read_label_file
        from Betsy import module_utils
        from genomicode import jmath
        import arrayio
        data_node, cls_node = antecedents
        # obtain the class label
        label, label_line, second_line = read_label_file.read(cls_node.identifier)
        class_num = len(label)
        assert class_num == 2, 'the number of class is not 2'
        fc = 1
        if 'group_fc_num' in user_options:
            fc = int(user_options['group_fc_num'])
        
        M = arrayio.read(data_node.identifier)
        first = M.slice(None, label[0][0])
        second = M.slice(None, label[1][0])
        #X = M.slice()
        I_good = []
        for i in range(M.nrow()):
            fold_change = abs(jmath.mean(first[i]) - jmath.mean(second[i]))
            if fold_change >= math.log(fc, 2):
                I_good.append(i)
        
        assert I_good, 'there is no gene is significant in fold change with 2'
        f = file(outfile, 'w')
        M_c = M.matrix(I_good, None)
        arrayio.tab_delimited_format.write(M_c, f)
        f.close()


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, cls_node = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        filename = 'signal_group_fc_' + original_file + '.tdf'
        return filename
