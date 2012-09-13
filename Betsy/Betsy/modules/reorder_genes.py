#reorder_genes.py
import hash_method
import gene_ranking
import module_utils
import os
import arrayio
from genomicode import arrayannot
def run(parameters, objects, pipeline):
    single_object = get_identifier(parameters, objects)
    outfile = get_outfile(parameters, objects, pipeline)
    #read the gene order list
    gene_list_file=module_utils.find_object(parameters,
                                objects, 'gene_list_file', 'contents')
    assert os.path.exists(gene_list_file.identifier), (
        'cannot find gene_list_file %s' % gene_list_file.identifier)
    gene_list = open(gene_list_file.identifier, 'r').read().split()
    M = arrayio.read(single_object.identifier)
    x = arrayannot.identify_all_platforms_of_matrix(M)
    id = x[0][0]
    platform = x[0][1]
    chip = arrayannot.identify_platform_of_annotations(gene_list)
    signal_file = single_object.identifier
    if platform == chip:
        tmpfile = single_object.identifier
    else:
        if parameters['platform'] in[chip, 'unknown_platform']:
            import subprocess
            import config
            Annot_path = config.ANNOTATE_MATRIX
            Annot_BIN = module_utils.which(Annot_path)
            assert Annot_BIN, 'cannot find the %s' % Annot_path
            signal_file = 'tmp'    
            command = ['python', Annot_BIN, '-f', single_object.identifier,
                       '-o', signal_file, "--platform", chip]
            process = subprocess.Popen(command, shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
            error_message = process.communicate()[1]
            if error_message:
                raise ValueError(error_message)
            assert module_utils.exists_nz(signal_file), 'the platform conversion fails'
            id = parameters['platform']
        elif parameters['platform'] == platform:
            infile = gene_list_file.identifier
            f=file(infile, 'rU')
            genes=f.readlines()
            f.close()
            gene_list = module_utils.convert_gene_list_platform(
                genes, platform)
    #read the tdf signal file
    M = arrayio.read(signal_file)
    original_list = M._row_names[id]
    #get the order index and write to the outout file
    indexlist = gene_ranking.find_sorted_index(original_list, gene_list)
    M_new = M.matrix(indexlist, None)
    f = open(outfile, 'w')
    arrayio.tab_delimited_format.write(M_new, f)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for reorder_genes fails' % outfile)
    new_objects = get_newobjects(parameters, objects, pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters, single_object, pipeline, outfile)
    return new_objects


def make_unique_hash(identifier, pipeline, parameters):
    return module_utils.make_unique_hash(
        identifier, pipeline, parameters)


def get_outfile(parameters, objects, pipeline):
    single_object = get_identifier(parameters, objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_reorder_' + original_file + '.tdf'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
                  parameters, objects, 'signal_file', 'contents,preprocess')
    assert os.path.exists(single_object.identifier), (
        'the input file %s for reorder_genes does not exist'
        %single_object.identifier)
    return single_object


def get_newobjects(parameters, objects, pipeline):
    outfile = get_outfile(parameters,objects, pipeline)
    single_object = get_identifier(parameters, objects)
    new_objects = module_utils.get_newobjects(
        outfile, 'signal_file', parameters, objects, single_object)
    return new_objects
