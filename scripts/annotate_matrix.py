#!/usr/bin/env python

from genomicode import arrayannot,jmath,Matrix
import arrayio
import os
import argparse

def main():
    
    parser = argparse.ArgumentParser(
        description = 'annotate the matrix')
    parser.add_argument('-o',dest = 'outpath',type = str,
                        help = 'specify the outpath',default=None)
    parser.add_argument('-f',dest = 'input',type = str,
                        help = 'specify the input file',default=None)
    parser.add_argument('--platform',dest = 'platform',type = str,
                        help = 'specify the platform to add',default=None)
    args = parser.parse_args()
    cwd = os.getcwd()
    platform = args.platform
    assert platform in platform2attributes,'we cannot convert to the platform %s'%platform
    in_id,in_platform = arrayannot.guess_platform(args.input)
    assert in_id, 'we cannot guess the platform for the input file'
    in_attribute,in_mart = platform2attributes[in_platform]
    out_attribute,out_mart = platform2attributes[platform]
    M = arrayio.read(args.input)
    gene_id_old = M._row_names[in_id]
    gene_id = ['"'+i+'"' for i in gene_id_old]
    R = jmath.start_R()
    jmath.R_equals_vector(gene_id,'gene_id')
    R('library(biomaRt)')
    jmath.R_equals('"'+in_attribute+'"','in_attribute')
    jmath.R_equals('"'+in_attribute+'"','filters')
    jmath.R_equals('"'+in_mart+'"','in_mart')
    R('old=useMart("ensembl",in_mart)')
    jmath.R_equals('"'+out_attribute+'"','out_attribute')
    jmath.R_equals('"'+out_mart+'"','out_mart')
    R('new=useMart("ensembl",out_mart)')
    R('homolog = getLDS(attributes=in_attribute,filters=filters,values=gene_id,mart=old,attributesL=out_attribute,martL=new)')
    homolog=R['homolog']
    old_id = [str(i) for i in homolog[0]]
    human_id = [str(i) for i in homolog[1]]
    new_id = []
    index_list=[]
    for i in range(len(gene_id_old)):
        gene_name = M._row_names[in_id][i]
        if gene_name in old_id:
            index_id = [index for index, item in enumerate(old_id) if item == gene_name]
            human_ids = [human_id[k] for k in index_id]
            new_id.extend(human_ids)
            index_old = gene_id_old.index(str(gene_name))
            index_list.extend([index_old]*len(index_id))
        else:
            new_id.append('')
            index_list.append(i)
    M_new = M.matrix(index_list,None)
    ids = M._row_order
    c=[platform]
    ids.extend(c)
    M_new._row_order=ids
    M_new._row_names[platform]=new_id
    newfile=file(args.outpath,'w')
    arrayio.tab_delimited_format.write(M_new,newfile)
    newfile.close()
    os.chdir(cwd)

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


    
if __name__=='__main__':
    main()
