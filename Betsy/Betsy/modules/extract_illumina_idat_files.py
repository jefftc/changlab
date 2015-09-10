from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        import shutil
        from genomicode import filelib
        from Betsy import module_utils
        
        path = module_utils.unzip_if_zip(in_data.identifier)
        x = filelib.list_files_in_path(path)
        x = [x for x in x if x.lower().endswith(".idat")]
        assert x, "No idat files."
        in_filenames = x

        if not os.path.exists(out_path):
            os.mkdir(out_path)
        for in_filename in in_filenames:
            in_path, in_file = os.path.split(in_filename)
            file_, ext = os.path.splitext(in_file)
            if file_.endswith("_Grn"):
                file_ = file_[:-4]
            out_file = "%s%s" % (file_, ext)
            out_filename = os.path.join(out_path, out_file)
            shutil.copyfile(in_filename, out_filename)
            
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
        return "idat_files"
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #filename = 'idat_files_' + original_file
        #return filename



