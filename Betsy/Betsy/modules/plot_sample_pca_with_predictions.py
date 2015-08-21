from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import genesetlib
        from Betsy import module_utils
        import plot_sample_pca
        data_node, classify_node = antecedents
        result_data = genesetlib.read_tdf(classify_node.identifier,
                                          preserve_spaces=True,
                                          allow_duplicates=True)
        for i in result_data:
            if i[0] == 'Predicted_class':
                legend = i[2]
        
        
        colors = ['r', 'b', 'g', 'y']
        legend_dict = {}
        for index, item in enumerate(legend):
            if item not in legend_dict:
                legend_dict[item] = [index]
            else:
                legend_dict[item].append(index)
        
        
        color = [''] * len(legend)
        for index, key in enumerate(legend_dict.keys()):
            c = colors[index]
            for i in legend_dict[key]:
                color[i] = c
        
        
        plot_sample_pca.plot_pca(data_node.identifier, outfile, color, legend)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, classify_node = antecedents
        original_file = module_utils.get_inputid(classify_node.identifier)
        loocv = ''
        if classify_node.data.attributes['loocv'] == 'yes':
            loocv = 'loocv'
        
        filename = ('prediction_pca_plot' + original_file + '_' +
                    classify_node.data.attributes['classify_alg'] + loocv + '.png')
        return filename


    def hash_input(self, pipeline, antecedents, out_attributes, user_options):
        from Betsy import module_utils
        data_node, classify_node = antecedents
        identifier = data_node.identifier
        return module_utils.hash_input(identifier, pipeline, out_attributes,
                                             user_options)


    def find_antecedents(
        self, network, module_id, out_attributes, user_attributes, pool):
        from Betsy import module_utils
        filter1 = module_utils.AntecedentFilter(datatype_name='PcaAnalysis')
        filter2 = module_utils.AntecedentFilter(
            datatype_name='ClassLabelFile',
            classify_alg=out_attributes["classify_alg"])
        x = module_utils.find_antecedents(
            network, module_id, user_attributes, pool, filter1, filter2)
        return x
