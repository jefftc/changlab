#gene_filter.py
import os
import module_utils

def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    import arrayio
    f_out = file(outfile,'w')
    M = arrayio.read(identifier)
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
    arrayio.pcl_format.write(M_c,f_out)
    f_out.close()
    assert module_utils.exists_nz(outfile),'the output file %s\
                                            for gene_filter fails'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'signal_file','contents',pipeline)
    
def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents')
    assert os.path.exists(identifier),'the input file %s\
                          for gene_filter does not exist'%identifier
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects
