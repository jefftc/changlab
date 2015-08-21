from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import filelib
        import os
        import arrayio
        from Betsy import module_utils
        from genomicode import config
        from genomicode import arrayplatformlib
        in_data = antecedents
        mapfile = config.HumanHT_12_to_HG_u133_Plus_2
        assert os.path.exists(mapfile), 'mapping file %s does not exist' % mapfile
        result = []
        for d in filelib.read_row(mapfile, header=True):
            if int(d.Distance) <= 1000 and d.Match == 'Best for Both':
                result.append((d.Affymetrix_Probe_Set_ID, d.Illumina_Probe_ID))
    
        
        
        M = arrayio.read(in_data.identifier)
        #platform_list = arrayplatformlib.identify_all_platforms_of_matrix(M)
        platform_list = arrayplatformlib.score_all_platforms_of_matrix(M)
        illu_id = None
        probe_id = None
        for platform in platform_list:
            if 'HumanHT_12' in platform:
                illu_id = M._row_names[platform[0]]
            if 'HG_U133_Plus_2' in platform:
                probe_id = M._row_names[platform[0]]
    
        
        
        if not illu_id or not probe_id:
            return None
    
        
        
        index = []
        for i in range(M.nrow()):
            if (probe_id[i], illu_id[i]) in result:
                index.append(i)
    
        
        
        if len(index) > 0:
            M_new = M.matrix(index, None)
            f = file(outfile, 'w')
            arrayio.tab_delimited_format.write(M_new, f)
            f.close()
            assert module_utils.exists_nz(outfile), (
                'the output file %s for best_match_both fails' % outfile
            )
        else:
            return None




    
    
    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'signal_best_match_both_' + original_file + '.tdf'
        return filename



