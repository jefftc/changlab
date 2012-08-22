#!/usr/bin/env python

from genomicode import arrayannot,jmath,Matrix,arrayplatformlib
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
    out_platform = args.platform
    assert arrayplatformlib.get_bm_organism(out_platform),'we cannot convert to the platform %s'%platform
    DATA = arrayio.read(args.input)
    platform_list = arrayannot.identify_all_platforms_of_matrix(DATA)
    assert platform_list, 'we cannot guess the platform for the input file'
    in_id = platform_list[0][0]
    in_platform = platform_list[0][1]
    in_attribute = arrayplatformlib.get_bm_attribute(in_platform)
    in_mart = arrayplatformlib.get_bm_organism(in_platform)
    out_attribute = arrayplatformlib.get_bm_attribute(out_platform)
    out_mart = arrayplatformlib.get_bm_organism(out_platform)
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
    c=[out_platform]
    ids.extend(c)
    M_new._row_order=ids
    M_new._row_names[out_platform]=new_id
    newfile=file(args.outpath,'w')
    arrayio.tab_delimited_format.write(M_new,newfile)
    newfile.close()
    os.chdir(cwd)

    
if __name__=='__main__':
    main()
