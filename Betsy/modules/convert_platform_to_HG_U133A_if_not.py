#convert_platform_to_affyu1332_if_not.py

import os
import module_utils
import subprocess
import Betsy_config
from genomicode import arrayannot,jmath,Matrix
import shutil
import arrayio

Databases = ["hsapiens_gene_ensembl","mmusculus_gene_ensembl"]

def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    slice_path = Betsy_config.RENAME
    slice_BIN = module_utils.which(slice_path)
    assert slice_BIN,'cannot find the %s' %slice_path
    chipname = arrayannot.guess_chip(single_object.identifier)
    if chipname == 'hg_u133_plus_2':
        shutil.copyfile(single_object.identifier,outfile)
    else:
        attributes,mart = check_platform_in_biomart(chipname)
        if attributes:
            M = arrayio.read(single_object.identifier)
            M_new = convert_others_to_hg_u133plus2(M,attributes,mart)
            f = file(outfile,'w')
            arrayio.tab_delimited_format.write(M_new,f)
            f.close()
        elif not attributes:
            raise ValueError('we cannot convert the platform you input to HG_U133plus2')
    assert module_utils.exists_nz(outfile),'the output file %s\
                                            for convert_platform fails'%outfile
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
    return new_objects


def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'signal' + '_u133plus2_'+original_file+'.pcl'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),'the input file %s\
                          for convert_platform'%single_object.identifier
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'signal_file',parameters,objects,single_object)
    return new_objects

def convert_others_to_hg_u133plus2(M,attributes,mart):
    ids = M._row_order
    gene_id_old = M._row_names[ids[0]]
    gene_id = ['"'+i+'"' for i in gene_id_old]
    R = jmath.start_R()
    jmath.R_equals_vector(gene_id,'gene_id')
    R('library(biomaRt)')
    R('human=useMart("ensembl","hsapiens_gene_ensembl")')
    command = 'old=useMart("ensembl",'+'"'+mart+'")'
    R(command)
    jmath.R_equals('"'+attributes+'"','filters')
    jmath.R_equals('"'+attributes+'"','attributes')
    R('homolog = getLDS(attributes=attributes,filters=filters,values=gene_id,mart=old,attributesL="affy_hg_u133_plus_2",martL=human)')
    homolog=R['homolog']
    old_id = [i for i in homolog[0]]
    human_id = [i for i in homolog[1]]
    index_list=[]
    for gene_name in old_id:
        index = gene_id_old.index(gene_name)
        index_list.append(index)
    M_new = M.matrix(index_list,None)
    ids = M_new._row_order
    c=['Probe.Set.ID']
    c.extend(ids)
    M_new._row_order=c
    M_new._row_names['Probe.Set.ID']=human_id
    return M_new

def check_platform_in_biomart(chipname):
    R = jmath.start_R()
    R('library(biomaRt)')
    for db in Databases:
        attributes = check_one_mart(db,chipname)
        if attributes:
            return attributes,db
    raise ValueError('we cannot find the platform')

def check_one_mart(ensembl,chipname):
    R = jmath.start_R()
    R('library(biomaRt)')
    command = 'db=useMart("ensembl",'+'"'+ensembl+'")'
    R(command)
    R('a=listAttributes(db)')
    a = R['a'][0]
    attributes = None
    # the match the chipname in mart e.g."affy_hg_u133a" with "hg_u133a"
    for j in a:
        jj = j.split('_')
        jj = '_'.join(jj[1:])
        if chipname.lower()==jj:
            attributes = j 
            mart = 'mouse'
    return attributes


