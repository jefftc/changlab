# TODO: Check if these are defined somewhere else.
MANIFEST_ALL = [
    'HumanHT-12_V3_0_R2_11283641_A.txt',
    'HumanHT-12_V4_0_R2_15002873_B.txt',
    'HumanHT-12_V3_0_R3_11283641_A.txt',
    'HumanHT-12_V4_0_R1_15002873_B.txt',
    'HumanMI_V1_R2_XS0000122-MAP.txt',
    'HumanMI_V2_R0_XS0000124-MAP.txt',
    'HumanRef-8_V2_0_R4_11223162_A.txt',
    'HumanRef-8_V3_0_R1_11282963_A_WGDASL.txt',
    'HumanRef-8_V3_0_R2_11282963_A.txt',
    'HumanRef-8_V3_0_R3_11282963_A.txt',
    'HumanWG-6_V2_0_R4_11223189_A.txt',
    'HumanWG-6_V3_0_R2_11282955_A.txt',
    'HumanWG-6_V3_0_R3_11282955_A.txt',
    'MouseMI_V1_R2_XS0000127-MAP.txt',
    'MouseMI_V2_R0_XS0000129-MAP.txt',
    'MouseRef-8_V1_1_R4_11234312_A.txt',
    'MouseRef-8_V2_0_R2_11278551_A.txt',
    'MouseRef-8_V2_0_R3_11278551_A.txt',
    'MouseWG-6_V1_1_R4_11234304_A.txt',
    'MouseWG-6_V2_0_R2_11278593_A.txt',
    'MouseWG-6_V2_0_R3_11278593_A.txt',
    'RatRef-12_V1_0_R5_11222119_A.txt',
    ]

CHIP_ALL = [
    'ilmn_HumanHT_12_V3_0_R3_11283641_A.chip',
    'ilmn_HumanHT_12_V4_0_R1_15002873_B.chip',
    'ilmn_HumanRef_8_V2_0_R4_11223162_A.chip',
    'ilmn_HumanReF_8_V3_0_R1_11282963_A_WGDASL.chip',
    'ilmn_HumanRef_8_V3_0_R3_11282963_A.chip',
    'ilmn_HumanWG_6_V2_0_R4_11223189_A.chip',
    'ilmn_HumanWG_6_V3_0_R3_11282955_A.chip',
    'ilmn_MouseRef_8_V1_1_R4_11234312_A.chip',
    'ilmn_MouseRef_8_V2_0_R3_11278551_A.chip',
    'ilmn_MouseWG_6_V1_1_R4_11234304_A.chip',
    'ilmn_MouseWG_6_V2_0_R3_11278593_A.chip',
    'ilmn_RatRef_12_V1_0_R5_11222119_A.chip',
    ]

from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import os
        import subprocess
        import arrayio
        from Betsy import module_utils
        from genomicode import config
        import zipfile
    
        in_data = antecedents
    
        module_name = 'IlluminaExpressionFileCreator'
        gp_parameters = dict()
        if zipfile.is_zipfile(in_data.identifier):
            gp_parameters['idat.zip'] = in_data.identifier
        else:
            # Add ".zip" to the end of the file.
            zipfile_name = os.path.split(in_data.identifier)[-1] + '.zip'
            zip_directory(in_data.identifier, zipfile_name)
            gp_parameters['idat.zip'] = os.path.join(".", zipfile_name)
    
        
        gp_parameters['manifest'] = 'HumanHT-12_V4_0_R2_15002873_B.txt'
        gp_parameters['chip'] = 'ilmn_HumanHT_12_V4_0_R1_15002873_B.chip'
        if 'illu_bg_mode' in out_attributes.keys():
            assert out_attributes['illu_bg_mode'] in [
                'true', 'false'
            ], 'illu_bg_mode should be ill_yes or ill_no'
            gp_parameters['background.subtraction.mode'] = \
                                                         out_attributes['illu_bg_mode']

        if 'illu_chip' in out_attributes.keys():
            assert out_attributes['illu_chip'] in CHIP_ALL, \
                   'illu_chip is not correct'
            gp_parameters['chip'] = str(out_attributes['illu_chip'])
        if 'illu_manifest' in out_attributes.keys():
            assert out_attributes['illu_manifest'] in MANIFEST_ALL, \
                   'illu_manifest is not correct'
            gp_parameters['manifest'] = str(out_attributes['illu_manifest'])
        if 'illu_coll_mode' in out_attributes.keys():
            assert out_attributes['illu_coll_mode'] in [
                'none', 'max', 'median'
            ], 'ill_coll_mode is not correct'
            gp_parameters['collapse.mode'] = str(
                out_attributes['illu_coll_mode'])
        if 'illu_clm' in user_options:
            gp_parameters['clm'] = str(user_options['illu_clm'])
        if 'illu_custom_chip' in user_options:
            gp_parameters['chip'] = str(user_options['illu_custom_chip'])
        if 'illu_custom_manifest' in user_options:
            gp_parameters['custom.manifest'] = str(
                user_options['illu_custom_manifest'])
    
        
        gp_path = config.genepattern
        gp_module = module_utils.which(gp_path)
        assert gp_module, 'cannot find the %s' % gp_path
        download_directory = os.path.join(".", 'illumina_result')
        command = [gp_module, module_name, '-o', download_directory]
        for key in gp_parameters.keys():
            a = ['--parameters', key + ':' + gp_parameters[key]]
            command.extend(a)
    
        
        process = subprocess.Popen(command,
                                   shell=False,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        process.wait()
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
    
        
        #goal_file = None
        assert os.path.exists(download_directory), (
            'there is no output directory for illumina'
        )
        result_files = os.listdir(download_directory)
        #assert 'stderr.txt' not in result_files,('gene_pattern get error '+
        #        'The contents of stderr.txt is:'+
        #        file(os.path.join(download_directory,'stderr.txt')).read())

        if not os.path.exists(outfile):
            os.mkdir(outfile)
    
        
        for result_file in result_files:
            if result_file == 'System.out':
                continue
            result_path = os.path.join(outfile, result_file)
            M = arrayio.read(os.path.join(download_directory, result_file))
            a = M._col_names['_SAMPLE_NAME']
            b = sorted(a)
            index = []
            for i in b:
                index.append(a.index(i))
            M_new = M.matrix(None, index)
            f = file(result_path, 'w')
            arrayio.gct_format.write(M_new, f)
            f.close()
    
        
        assert module_utils.exists_nz(outfile), (
            'the output file %s for illumina fails' % outfile)
        

    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'IlluFolder_' + original_file
        return filename


def zip_directory(dir, zip_file):
    import os
    import zipfile

    zip = zipfile.ZipFile(zip_file, 'w', compression=zipfile.ZIP_DEFLATED)
    root_len = len(os.path.abspath(dir))
    for root, dirs, files in os.walk(dir):
        archive_root = os.path.abspath(root)[root_len:]
        for f in files:
            fullpath = os.path.join(root, f)
            archive_name = os.path.join(archive_root, f)
            zip.write(fullpath, archive_name, zipfile.ZIP_DEFLATED)
    
    zip.close()
