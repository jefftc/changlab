from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import os
        import subprocess
        import shutil
        import arrayio
        from Betsy import module_utils
        from genomicode import config
        from genomicode import arrayplatformlib
        data_node, cls_node = antecedents
        gsea_path = config.gsea
        gsea_module = module_utils.which(gsea_path)
        assert gsea_module, 'cannot find the %s' % gsea_path
        M = arrayio.read(data_node.identifier)
        x = arrayplatformlib.identify_all_platforms_of_matrix(M)
        chipname = x[0][1]
        assert chipname in arrayplatformlib.platform_to_GSEA_chipname, (
            'we cannot find chipname %s in gsea' % chipname
        )
        chipname_in_gsea = arrayplatformlib.platform_to_GSEA_chipname[chipname]
        platform = chipname_in_gsea + '.chip'
        download_directory = os.path.join(".", 'gsea_result')
        command = [gsea_module, data_node.identifier, '--cls_file',
                   cls_node.identifier, '--platform', platform, '--database',
                   out_attributes['geneset_database'], '--permutation_type',
                   out_attributes['permutation_type'], download_directory]
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
        
        assert module_utils.exists_nz(download_directory), (
            'there is no output directory for GSEA'
        )
        shutil.copytree(download_directory, outfile)
        assert module_utils.exists_nz(outfile), (
            'the output file %s for annotate_genes_with_gsea fails' % outfile
        )


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node, cls_node = antecedents
        original_file = module_utils.get_inputid(data_node.identifier)
        filename = 'gsea_' + original_file
        return filename


    
