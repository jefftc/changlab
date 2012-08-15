#add_HG_U133A_probeids.py

import os
import module_utils
import subprocess
import Betsy_config
from genomicode import arrayannot,jmath,Matrix
import shutil
import arrayio

def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    chipname = arrayannot.guess_chip(single_object.identifier)
    if chipname == 'HG_U133_Plus_2':
        shutil.copyfile(single_object.identifier,outfile)
    else:
        if chipname in platform2attributes:
            attributes,mart = platform2attributes[chipname]
            M = arrayio.read(single_object.identifier)
            M_new = convert_others_to_hg_u133plus2(M,attributes,mart)
            f = file(outfile,'w')
            arrayio.tab_delimited_format.write(M_new,f)
            f.close()
        else:
            raise ValueError('we cannot convert the platform you input to HG_U133plus2')
    assert module_utils.exists_nz(outfile),(
        'the output file %s for convert_platform fails'%outfile)
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
    assert os.path.exists(single_object.identifier),(
        'the input file %s for convert_platform'%single_object.identifier)
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
    #R('homolog = getLDS(attributes=attributes,filters=filters,values=gene_id,mart=old,attributesL="entrezgene",martL=human)')
    homolog=R['homolog']
    old_id = [i for i in homolog[0]]
    human_id = [i for i in homolog[1]]
    index_list=[]
    for gene_name in old_id:
        index = gene_id_old.index(str(gene_name))
        index_list.append(index)
    M_new = M.matrix(index_list,None)
    c=['Affymetrix.Probe.Set.ID']
    c.extend(ids)
    M_new._row_order=c
    M_new._row_names['Affymetrix.Probe.Set.ID']=human_id
    return M_new

platform2attributes={
                 'HG_U133_Plus_2':("affy_hg_u133_plus_2","hsapiens_gene_ensembl"),
                 'HG_U133B':("affy_hg_u133b","hsapiens_gene_ensembl"),
                 'HG_U133A':("affy_hg_u133a","hsapiens_gene_ensembl"),
                 'HG_U133A_2':("affy_hg_u133a_2","hsapiens_gene_ensembl"),
                 'HG_U95A':("affy_hg_u95a","hsapiens_gene_ensembl"),
                 'HumanHT_12':("illumina_humanht_12","hsapiens_gene_ensembl"),
                 'HG_U95Av2':("affy_hg_u95av2","hsapiens_gene_ensembl"),
                 'entrez_ID_human':("entrezgene","hsapiens_gene_ensembl"),
                 'entrez_ID_symbol_human':("hgnc_symbol","hsapiens_gene_ensembl"),
                 'Hu6800':("affy_hugenefl","hsapiens_gene_ensembl"),
                 
                 'Mouse430A_2':('affy_mouse430a_2',"mmusculus_gene_ensembl"),
                 'MG_U74Cv2':('affy_mg_u74cv2',"mmusculus_gene_ensembl"),
                 'Mu11KsubB':("affy_mu11ksubb","mmusculus_gene_ensembl"),
                 'Mu11KsubA':('affy_mu11ksuba',"mmusculus_gene_ensembl"),
                 'MG_U74Av2':("affy_mg_u74av2","mmusculus_gene_ensembl"),
                 'Mouse430_2':('affy_mouse430_2',"mmusculus_gene_ensembl"),
                 'MG_U74Bv2':('affy_mg_u74bv2',"mmusculus_gene_ensembl"),
                 'entrez_ID_mouse':("entrezgene","mmusculus_gene_ensembl"),
                 'MouseRef_8':("illumina_mousewg_6_v2","mmusculus_gene_ensembl"),
                 'entrez_ID_symbol_mouse':("mgi_symbol","mmusculus_gene_ensembl"),
                 
                 'RG_U34A':('affy_rg_u34a',"rnorvegicus_gene_ensembl"),
                 'RAE230A':('affy_rae230a',"rnorvegicus_gene_ensembl")}

                

