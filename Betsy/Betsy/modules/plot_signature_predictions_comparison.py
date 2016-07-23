from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """generate a heatmap of input file"""
        from genomicode import graphlib
        from genomicode import filelib
        import subprocess
        import arrayio
        from Betsy import module_utils
        from genomicode import config
        in_data = antecedents
        Heatmap_path = config.arrayplot
        Heatmap_BIN = module_utils.which(Heatmap_path)
        assert Heatmap_BIN, 'cannot find the %s' % Heatmap_path

        command = ['python', Heatmap_BIN, in_data.identifier, '-o', outfile,
                   "--label_arrays", "--label_genes"]
        if 'color' in out_attributes.keys():
            color = ['--color', out_attributes['color'].replace('_', '-')]
            command.extend(color)
    
        
        
        M = arrayio.read(in_data.identifier)
        nrow = M.nrow()
        ncol = M.ncol()
        ratio = float(nrow) / ncol
        max_box_height = 20
        max_box_width = 60
        if 'hm_width' in user_options:
            max_box_width = user_options['hm_width']
    
        
        
        if 'hm_height' in user_options:
            max_box_height = user_options['hm_height']
    
        
        
        if ratio >= 4:
            x, y = graphlib.find_tall_heatmap_size(nrow, ncol,
                                                   max_box_height=max_box_height,
                                                   max_box_width=max_box_width,
                                                   min_box_height=20,
                                                   min_box_width=20,
                                                   max_megapixels=128)
        else:
            x, y = graphlib.find_wide_heatmap_size(nrow, ncol,
                                                   max_box_height=max_box_height,
                                                   max_box_width=max_box_width,
                                                   min_box_height=20,
                                                   min_box_width=20,
                                                   max_megapixels=128)
    
        
        
        command.extend(['-x', str(x), '-y', str(y)])
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
    
        
        
        assert filelib.exists_nz(outfile), (
            'the output file %s for plot_signature_prediction_comparison fails' %
            outfile)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'heatmap_' + original_file + '.png'
        return filename



