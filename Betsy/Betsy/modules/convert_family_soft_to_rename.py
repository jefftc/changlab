from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        """convert family soft file to RenameFile"""
        from Betsy import module_utils
        in_data = antecedents
        GSEID = user_options['GSEID']
        title, description = extract_sample2desc(GSEID, in_data.identifier)
        if out_attributes['labels_from'] == 'title':
            convert_family_to_relabel_file(title, outfile)
        else:
            convert_family_to_relabel_file(description, outfile)
    
        assert module_utils.exists_nz(outfile), (
            'the output file %s for convert_family_soft_to_rename does not exists'
            % outfile)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'rename_file_' + original_file + '.tdf'
        return filename


def extract_sample2desc(GSEID, filename):
    from genomicode.filelib import openfh
    title_dict = {}
    description_dict = {}
    id_ = None
    for line in openfh(filename):
        if line.startswith("^SAMPLE"):
            assert id_ is None, "problem with %s" % filename
            id_ = line.strip().split()[2]
        elif line.startswith("!Sample_description"):
            assert id_ is not None, "problem with %s" % filename
            title = line.strip().split(None, 2)[2]
            #x = id, title
            description_dict[id_] = title
        elif line.startswith("!Sample_title"):
            assert id_ is not None, "problem with %s" % filename
            title = line.strip().split(None, 2)[2]
            #x = id, "Title: %s" % title
            title_dict[id_] = title
        elif line.startswith("!sample_table_end"):
            id_ = None
    
    return title_dict, description_dict


def convert_family_to_relabel_file(relabel_dict, outfile):
    f = file(outfile, 'w')
    f.write('\t'.join(['Name', 'NewName']) + '\n')
    for key in relabel_dict:
        f.write('\t'.join([key, relabel_dict[key]]) + '\n')
    f.close()
    
