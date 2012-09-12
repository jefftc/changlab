#filter_genes_by_missing_values.py
import os
import module_utils

def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    import arrayio
    f_out = file(outfile,'w')
    M = arrayio.read(single_object.identifier)
    I_good = []
    #get the percentage of gene filter
    percent=float(parameters['filter'])/100
    for i in range(M.dim()[0]):
       missing_count = 0
       for j in range(M.dim()[1]):
           if M._X[i][j] == None:
                missing_count = missing_count + 1
       if float(missing_count)/M.dim()[1]<percent:
            I_good.append(i)
    M_c = M.matrix(I_good,None)
    arrayio.tab_delimited_format.write(M_c,f_out)
    f_out.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for gene_filter fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal_filter_'+original_file+'.tdf'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for gene_filter does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
