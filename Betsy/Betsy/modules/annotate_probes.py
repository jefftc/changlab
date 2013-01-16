#annotate_probes.py
from genomicode import arrayannot,arrayplatformlib
import arrayio
from Betsy import module_utils
import os

def run(parameters, objects, pipeline):
    single_object = get_identifier(parameters, objects)
    outfile = get_outfile(parameters, objects, pipeline)
    M = arrayio.read(single_object.identifier)
    all_platforms = arrayplatformlib.identify_all_platforms_of_matrix(M)
    if not all_platforms:
        raise ValueError('we cannot guess the platform and annotate the file')
    ids = M._row_order
    probe_header = all_platforms[0][0]
    probe_id = M._row_names[probe_header]
    new_ids = ids[:]
    new_ids.remove(probe_header)
    annotate_type = parameters['annotate_type']
    if annotate_type == 'all':
        annotate_header = arrayplatformlib.annotate_header
    elif annotate_type == 'gene_id':
        annotate_header = ['Gene ID']
    dictionary = arrayannot.annotate_probes_multiple(probe_id, annotate_header)
    column=[]
    for id in new_ids:
        flag = True
        for key in dictionary.keys():
            flag = True
            value_list = dictionary[key]
            gene_list = M._row_names[id]
            for gene in gene_list:
                if gene not in value_list:
                    flag = False
                    break
            if flag:
                column.append((key,id))
    header = [i[0] for i in column]
    miss_header = list(set(annotate_header).difference(set(header)))
    original_ids = ids[:]
    for col in miss_header:
        col_2 = col
        if col in original_ids:
            col_1 = col + '_1'
            col_2 = col + '_2'
            M, ids = module_utils.replace_matrix_header(M,col,col_1)
        ids.append(col_2)
        M._row_order = ids
        M._row_names[col_2] = dictionary[col]
    f = file(outfile, 'w')
    arrayio.tab_delimited_format.write(M, f)
    f.close()
    assert module_utils.exists_nz(outfile), (
        'the output file %s for annot_probes fails' % outfile)
    new_objects = get_newobjects(parameters, objects, pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters, single_object, pipeline, outfile)
    return new_objects


def make_unique_hash(identifier, pipeline, parameters):
    return module_utils.make_unique_hash(identifier, pipeline, parameters)


def get_outfile(parameters, objects, pipeline):
    single_object = get_identifier(parameters, objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = parameters['filetype'] + '_annot_' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_identifier(parameters, objects):
    single_object = module_utils.find_object(
        parameters, objects, parameters['filetype'], 'contents,preprocess')
    assert os.path.exists(single_object.identifier), (
        'the input file %s for annot_probes does not exist'
        % single_object.identifier)
    return single_object


def get_newobjects(parameters, objects, pipeline):
    outfile = get_outfile(parameters, objects, pipeline)
    single_object = get_identifier(parameters, objects)
    new_objects = module_utils.get_newobjects(
        outfile, parameters['filetype'], parameters, objects, single_object)
    return new_objects
