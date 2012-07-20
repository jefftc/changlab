#convert_platform_to_affyu1332_if_not.py

import os
import module_utils
import subprocess
import Betsy_config
<<<<<<< HEAD
from genomicode import arrayannot,jmath,Matrix
import shutil
import arrayio
=======
from genomicode import arrayannot
import shutil
>>>>>>> a28b8e4a09c20c623e5321ef81cffaccaf7b15f8
def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    slice_path = Betsy_config.RENAME
    slice_BIN = module_utils.which(slice_path)
    assert slice_BIN,'cannot find the %s' %slice_path
    chipname1 = arrayannot.guess_chip(single_object.identifier)
    attributes,mart = check_platform_in_biomart(chipname)
    if attributes:
        M = arrayio.read(single_object.identifier)
        M_new = convert_mouse_to_human(M,attributes,mart)
        f = file(outfile,'w')
        arrayio.tab_delimited_format.write(M_new,f)
        f.close()
<<<<<<< HEAD
    elif not attributes:
        if chipname1 in ['HumanHT-12']:
            mapfile = Betsy_config.MAPPING
            assert os.path.exists(mapfile),'the mapping file %s does not exists'%mapfile
            
            M = arrayio.read(single_object.identifier)
            num_header = len(M._row_names.keys())
            m = str(num_header+1)
            command1 = ['python', slice_BIN, '--select_row_annotation',mapfile+',Match,Best for Both',
                        '--select_row_numeric_annotation',mapfile+',Distance,<=1000',
                        '--add_row_annot',mapfile+',Affymetrix.Probe.Set.ID',
                         '--rename_row_annot', 'Affymetrix.Probe.Set.ID,Probe.Set.ID','--move_row_annot',m+',1',
                        single_object.identifier]
            f=file(outfile,'w')
            process = subprocess.Popen(command1,shell=False,
                                    stdout=f,
                                    stderr=subprocess.PIPE)
            error_message = process.communicate()[1]
            if error_message:
                    raise ValueError(error_message)
            f.close()
        
        elif parameters['preprocess'] in ['mas5','rma']:
                shutil.copyfile(single_object.identifier,outfile)
        else:
            raise ValueError('we cannot convert the platform you input to HG_U133A')
=======
    elif parameters['preprocess'] in ['mas5','rma']:
        shutil.copyfile(single_object.identifier,outfile)
    else:
        raise ValueError('we cannot convert the platform you input to HG_U133A')
>>>>>>> a28b8e4a09c20c623e5321ef81cffaccaf7b15f8
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
    filename = 'signal' + '_affyu1332_'+original_file+'.pcl'
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

def convert_mouse_to_human(M,attributes,mart):
    ids = M._row_order
    gene_id_old = M._row_names[ids[0]]
    gene_id = ['"'+i+'"' for i in gene_id_old]
    R = jmath.start_R()
    jmath.R_equals_vector(gene_id,'gene_id')
    R('library(biomaRt)')
    R('human=useMart("ensembl","hsapiens_gene_ensembl")')
    R('mouse=useMart("ensembl","mmusculus_gene_ensembl")')
    jmath.R_equals(attributes,'filters')
    jmath.R_equals(attributes,'attributes')
    jmath.R_equals(mart,'mart')
    R('homolog = getLDS(attributes=attributes,filters=filters,values=gene_id,mart=mart,attributesL="affy_hg_u133_plus_2",martL=human)')
    homolog=R['homolog']
    mouse_id = [i for i in homolog[0]]
    human_id = [i for i in homolog[1]]
    mouse_id = [mouse_id[i] for i in range(len(human_id)) if len(human_id[i])>0]
    human_id = [human_id[i] for i in range(len(human_id)) if len(human_id[i])>0]
    index_list=[]
    for gene_name in mouse_id:
        index = gene_id_old.index(gene_name)
        index_list.append(index)
    M_new = M.matrix(index_list,None)
    ids = M_new._row_order
    c=['gene_id']
    c.extend(ids)
    M_new._row_order=c
    M_new._row_names['gene_id']=human_id
    return M_new

def check_platform_in_biomart(chipname):
    R = jmath.start_R()
    R('library(biomaRt)')
    R('human=useMart("ensembl","hsapiens_gene_ensembl")')
    R('mouse=useMart("ensembl","mmusculus_gene_ensembl")')
    R('a=listAttributes(mouse)')
    R('b=listAttributes(human)')
    a = R['a'][0]
    b = R['b'][0]
    mart = None
    attributes = None
    for j in a:
        jj = j.split('_')
        jj = '_'.join(jj[1:])
        if chipname.lower()==jj:
            attributes = '"'+j +'"'
            mart = 'mouse' 
    for k in b:
        kk = k.split('_')
        kk = '_'.join(kk[1:])
        if chipname.lower() == kk:
            attributes = '"'+k +'"'
            mart = 'human' 
    return attributes,mart

