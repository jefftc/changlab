from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """check bamfiles in bam folder is all sorted"""
        import shutil
        from Betsy import module_utils
        in_data = antecedents
        #filenames = os.listdir(in_data.identifier)
        #out_attributes = set_out_attributes(in_data, out_attributes)
        shutil.copytree(in_data.identifier, outfile)
        assert module_utils.exists_nz(outfile), (
            'the output file %s for is_Bam_folder_sorted fails' % outfile
        )



    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'BamFolder_' + original_file
        return filename


    def set_out_attributes(self, antecedents, out_attributes):
        import os
        new_parameters = out_attributes.copy()
        filenames = os.listdir(antecedents.identifier)
        assert filenames, 'The input folder is empty.'
        sorted_list = []
        for filename in filenames:
            if filename == '.DS_Store':
                pass
            else:
                fileloc = os.path.join(antecedents.identifier, filename)
                sort_info = get_sort_information(fileloc)
                sorted_list.append(sort_info)
        
        if False in sorted_list:
            new_parameters['sorted'] = 'no'
        else:
            new_parameters['sorted'] = 'yes'
        
        return new_parameters



def get_sort_information(filename):
    import os
    from genomicode import config
    import subprocess
    samtools_BIN = config.samtools
    assert os.path.exists(samtools_BIN), 'cannot find the %s' % samtools_BIN
    command = [samtools_BIN, 'view', '-H', filename]
    process = subprocess.Popen(command,
                               shell=False,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    process.wait()
    out_message = process.communicate()
    if 'error' in out_message[1]:
        error_message = "\n".join(out_message)
        raise ValueError(error_message)
    
    message = out_message[0]
    header = None
    for line in message.split('\n'):
        if line.startswith('@HD'):
            header = line
            break
    
    sort = header.split('SO:')[-1]
    if sort == 'coordinate':
        return True
    elif sort == 'unsorted':
        return False
    else:
        raise ValueError('cannot determine sorted or not')

