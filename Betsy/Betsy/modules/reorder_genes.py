from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import arrayio
        from genomicode import config
        from genomicode import arrayplatformlib
        from Betsy import gene_ranking
        from Betsy import module_utils
        
        data_node, gene_node = antecedents
        M = arrayio.read(data_node.identifier)
        gene_list = open(gene_node.identifier, 'r').read().split()

        I_row, I_col = M._index(row=gene_list)
        M_new = M.matrix(I_row, None)

##         scores = arrayplatformlib.score_matrix(M, min_score=0.80)
##         x = arrayplatformlib.identify_all_platforms_of_matrix(M)
##         if x:
##             id = x[0][0]
##             platform = x[0][1]
##             chip = arrayplatformlib.identify_platform_of_annotations(gene_list)
##             if not chip:
##                 chip = []
##             signal_file = data_node.identifier
##             #if platform == chip:
##             #    tmpfile = data_node.identifier
##             if platform != chip:
##                 platform_name = 'unknown_platform'
##                 if 'platform_name' in user_options:
##                     platform_name = user_options['platform_name']
##                 if platform_name in chip:  #, 'unknown_platform':
##                     import subprocess
##                     Annot_path = config.annotate_matrix
##                     Annot_BIN = module_utils.which(Annot_path)
##                     assert Annot_BIN, 'cannot find the %s' % Annot_path
##                     signal_file = 'tmp'
##                     command = [
##                         'python',
##                         Annot_BIN,
##                         # What is single_object?  Bug?
##                         #'-f', single_object.identifier,
##                         '-f', data_node.identifier,
##                         '-o', signal_file,
##                         "--platform", chip,
##                         ]
##                     process = subprocess.Popen(command,
##                                                shell=False,
##                                                stdout=subprocess.PIPE,
##                                                stderr=subprocess.PIPE)
##                     error_message = process.communicate()[1]
##                     if error_message:
##                         raise ValueError(error_message)
##                     assert module_utils.exists_nz(
##                         signal_file), 'the platform conversion fails'
##                     id = out_attributes['platform']
##                     M = arrayio.read(signal_file)
##                 elif platform_name == platform:
##                     # Needs to be tested.
##                     #infile = gene_list_file.identifier
##                     infile = gene_node.identifier
##                     f = file(infile, 'rU')
##                     genes = f.readlines()
##                     f.close()
##                     gene_list = module_utils.convert_gene_list_platform(genes,
##                                                                         platform)
##         else:
##             id = M._row_order[0]
##         #id = M._row_order[3]
##         #print id
        
        
##         original_list = M._row_names[id]
##         #get the order index and write to the outout file
##         indexlist = gene_ranking.find_sorted_index(original_list, gene_list)
##         M_new = M.matrix(indexlist, None)
        
        f = open(outfile, 'w')
        arrayio.tab_delimited_format.write(M_new, f)
        f.close()


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, cls_node = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        filename = 'gene_reorder' + original_file + '.tdf'
        return filename


    
