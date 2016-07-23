from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import mplgraph
        import subprocess
        import arrayio
        from Betsy import module_utils
        from genomicode import config
        from genomicode import arrayplatformlib
        from genomicode import filelib
        
        in_data = antecedents
        M = arrayio.read(in_data.identifier)
        header = M._row_names
        header_name = M._row_order
        column_id = None
        for i in header_name:
            column = [j.upper() for j in header[i]]
            if 'ACTB' in column:
                column_id = i
                break
    
        
        if not column_id:
            x = arrayplatformlib.identify_all_platforms_of_matrix(M)
            #id = x[0][0]
            platform = x[0][1]
            if platform in [
                'HG_U133_Plus_2', 'HG_U133B', 'HG_U133A', 'HG_U133A_2',
                'HG_U95A', 'HumanHT_12', 'HG_U95Av2',
                'Entrez_ID_human', 'Entrez_Symbol_human', 'Hu6800']:
                out_platform = 'Entrez_Symbol_human'
            elif platform in [
                'Mouse430A_2', 'MG_U74Cv2', 'Mu11KsubB', 'Mu11KsubA',
                'MG_U74Av2', 'Mouse430_2', 'MG_U74Bv2',
                'Entrez_ID_mouse', 'MouseRef_8',
                'Entrez_Symbol_mouse']:
                out_platform = 'Entrez_Symbol_mouse'
            Annot_path = config.annotate_matrix
            Annot_BIN = module_utils.which(Annot_path)
            assert Annot_BIN, 'cannot find the %s' % Annot_path

            command = [
                'python',
                Annot_BIN,
                in_data.identifier, '--platform',
                out_platform,
                ]
            f = file('tmp', 'w')
            try:
                process = subprocess.Popen(
                    command, shell=False, stdout=f, stderr=subprocess.PIPE)
            finally:
                f.close()
            error_message = process.communicate()[1]
            if error_message:
                raise ValueError(error_message)
            assert filelib.exists_nz('tmp'), 'failed platform conversion'
            column_id = out_platform
            M = arrayio.read('tmp')
        
        label = M._col_names['_SAMPLE_NAME']
        keywords = ['ACTB', 'TUBB']
        lines = []
        data = []
        legend_name = []
        # Don't like this algorithm.  High chance of mis-labeling.
        for keyword in keywords:
            for i in range(M.dim()[0]):
                if M._row_names[column_id][i].upper() == keyword:
                    data.append(M.slice()[i])
                    x = "%s (%s)" % (keyword, M._row_names[header_name[0]][i])
                    legend_name.append(x)
        
        for i in range(len(data)):
            line = [(j, data[i][j]) for j in range(len(data[i]))]
            lines.append(line)
            fig = mplgraph.lineplot(
                *lines, legend=legend_name, box_label=label,
                ylim_min=0, ylabel='Gene Expression Value')
            fig.savefig(outfile)
    
        if not lines:
            import matplotlib.pyplot as plt
            plt.clf()
            plt.plot([0, 0, 0, 0])
            plt.title('no ACTB or TUBB probes are found')
            plt.savefig(outfile)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'actb_plot' + original_file + '.png'
        return filename


