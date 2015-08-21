from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import os
        import shutil
        from Betsy import module_utils
        in_data = antecedents
        directory = module_utils.unzip_if_zip(in_data.identifier)
        illumina_file = []
        filenames = os.listdir(directory)
        assert filenames, 'The input folder or zip file is empty.'
        for dirpath, dirnames, filenames in os.walk(directory):
            for filename in [f for f in filenames if f.endswith(".idat")]:
                illumina_file.append(os.path.join(dirpath, filename))
    
        
        if illumina_file:
            os.mkdir(outfile)
            for filename in illumina_file:
                newfilename = os.path.split(filename)[-1]
                if newfilename[:-5].endswith('_Grn'):
                    newfilename = newfilename[:-9] + newfilename[-5:]
                new_file = os.path.join(outfile, newfilename)
                shutil.copyfile(filename, new_file)
            assert module_utils.exists_nz(outfile), (
                'the output file %s for extract_illumina_idat_files fails' % outfile
            )
        else:
            assert ValueError('There is no illumina idat file in the input.')

            ##def run(data_node, parameters, user_input,network):
            ##    outfile = name_outfile(data_node,user_input)
            ##    directory = module_utils.unzip_if_zip(data_node.identifier)
            ##    illumina_file = []
            ##    filenames = os.listdir(directory)
            ##    assert filenames, 'The input folder or zip file is empty.'
            ##    for filename in filenames:
            ##        if filename in ['.DS_Store', '._.DS_Store', '.Rapp.history']:
            ##            continue
            ##        if filename.endswith('.idat'):
            ##            illumina_file.append(filename)
            ##    if illumina_file:
            ##        os.mkdir(outfile)
            ##        for filename in illumina_file:
            ##            if filename[:-5].endswith('_Grn'):
            ##                newfilename = filename[:-9] + filename[-5:]
            ##            else:
            ##                newfilename = filename
            ##            old_file = os.path.join(directory, filename)
            ##            new_file = os.path.join(outfile, newfilename)
            ##            shutil.copyfile(old_file, new_file)
            ##        assert module_utils.exists_nz(outfile), (
            ##            'the output file %s for extract_illumina_idat_files fails'
            ##            % outfile)
            ##        out_node = bie3.Data(rulebase.IDATFiles,**parameters)
            ##    	out_object = module_utils.DataObject(out_node,outfile)
            ##    	return out_object
            ##    else:
            ##        print 'There is no illumina idat file in the input.'
            ##        return None



    
    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'idat_files_' + original_file
        return filename



