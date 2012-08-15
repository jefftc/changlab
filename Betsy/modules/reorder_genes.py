#reorder_genes.py
import hash_method
import gene_ranking
import module_utils
import os
from genomicode import arrayannot
def run(parameters,objects,pipeline):
    #also if input is other kind of file
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    #read the gene order list
    gene_list_file=module_utils.find_object(parameters,
                                objects,'gene_list_file','contents')
    assert os.path.exists(gene_list_file.identifier),(
        'cannot find gene_list_file %s'%gene_list_file.identifier)
    gene_list = open(gene_list_file.identifier,'r').read().split()
    id,platform=arrayannot.guess_platform(single_object.identifier)
    chip = arrayannot.guess_chip_from_probesets(gene_list)
    if platfrom == chip:
        tmpfile = single_object.identifier
    else:
        import subprocess
        import Betsy_config
        Annot_path = Betsy_config.ANNOTATE_MATRIX
        Annot_BIN = module_utils.which(Annot_path)
        assert Annot_BIN,'cannot find the %s' %Annot_path
        if parameters['platform'] == chip:
            infile = single_object.identifier
            tmpfile = 'tmp'
            pf = chip
            command = ['python', Annot_BIN,'-f',single_object.identifier,'-o','tmp',"--platform",
               chip]
            process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
            error_message = process.communicate()[1]
            if error_message:
                raise ValueError(error_message)
            assert module_utils.exists_nz('tmp'),'the platform conversion fails'
        elif parameters['platform'] == platform:
            infile = gene_list_file.identifier
            tmpfile = 'tmp'
            pf = platform
        
        
        
    #read the pcl signal file
    f_signal= open(tmpfile,'r')
    content = f_signal.readlines()
    f_signal.close()
    #get the original gene list
    original_list = []
    for i in range(1,len(content)):
        original_list.append(content[i].split()[0])
    #get the order index and write to the outout file
    indexlist = gene_ranking.find_sorted_index(original_list,gene_list)
    f = open(outfile,'w')
    f.write(content[0])
    for i in range(len(indexlist)):
        f.write(content[indexlist[i]+1])
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for reorder_genes fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline,outfile)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_reorder_' + original_file + '.pcl'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
                  parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for reorder_genes does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
