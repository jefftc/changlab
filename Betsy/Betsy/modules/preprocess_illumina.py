from Betsy.rules import ArrayPlatforms as AP

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
        params = {}
        if zipfile.is_zipfile(in_data.identifier):
            params['idat.zip'] = in_data.identifier
        else:
            # Add ".zip" to the end of the file.
            zipfile_name = os.path.split(in_data.identifier)[-1] + '.zip'
            zip_directory(in_data.identifier, zipfile_name)
            params['idat.zip'] = os.path.join(".", zipfile_name)

        x = user_options.get("illu_manifest", AP.DEFAULT_MANIFEST)
        assert x in AP.ILLU_MANIFEST, "Unknown manifest: %s" % x
        params['manifest'] = x

        x = user_options.get("illu_chip", AP.DEFAULT_CHIP)
        assert x in AP.ILLU_CHIP, "Unknown chip: %s" % x
        params['chip'] = x

        x = user_options.get("illu_bg_mode")
        if x is not None:
            assert x in ['true', 'false'], \
                   'illu_bg_mode should be true or false'
            params['background.subtraction.mode'] = x

        x = user_options.get("illu_coll_mode")
        if x is not None:
            assert x in ['none', 'max', 'median'], \
                   'ill_coll_mode is not correct'
            params['collapse.mode'] = str(x)
            
        if 'illu_clm' in user_options:
            params['clm'] = str(user_options['illu_clm'])
        if 'illu_custom_chip' in user_options:
            params['chip'] = str(user_options['illu_custom_chip'])
        if 'illu_custom_manifest' in user_options:
            params['custom.manifest'] = str(
                user_options['illu_custom_manifest'])

        gp_path = config.genepattern
        gp_module = module_utils.which(gp_path)
        assert gp_module, 'cannot find: %s' % gp_path
        download_directory = 'illumina_result'
        command = [
            gp_module,
            module_name,
            '-o', download_directory,
            ]
        for key in params.keys():
            x = ['--parameters', key + ':' + params[key]]
            command.extend(x)


        process = subprocess.Popen(
            command, shell=False, stdout=subprocess.PIPE,
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
            # BUG: What if there are duplicate sample names?
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


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'IlluFolder_' + original_file
        return filename


# TODO: Move this to a library.
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


