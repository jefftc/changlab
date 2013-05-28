#!/usr/bin/env python

from genomicode import jmath,Matrix,arrayplatformlib
import arrayio
import os
import argparse

def main():
    all_platforms = [platform.name for platform in arrayplatformlib.PLATFORMS]
    parser = argparse.ArgumentParser(
        description = 'annotate the matrix')
    parser.add_argument('-o',dest = 'outpath',type = str,
                        help = 'specify the outpath',default=None)
    parser.add_argument('-f',dest = 'input',type = str,
                        help = 'specify the input file',default=None)
    parser.add_argument('--platform',dest = 'platform',type = str,
                        help = 'specify the platform to add:'+ str(all_platforms),
                        default=[],action = 'append')
    args = parser.parse_args()
    cwd = os.getcwd()
    out_platforms = args.platform
    assert out_platforms,'please give at least one platform for convert'
    for out_platform in out_platforms:
        assert arrayplatformlib.get_bm_organism(out_platform),'we cannot convert to the platform %s'%out_platform
    DATA = arrayio.read(args.input)
    platform_list = arrayplatformlib.identify_all_platforms_of_matrix(DATA)
    assert platform_list, 'we cannot guess the platform for the input file'
    in_id = platform_list[0][0]
    in_platform = platform_list[0][1]
    in_attribute = arrayplatformlib.get_bm_attribute(in_platform)
    in_mart = arrayplatformlib.get_bm_organism(in_platform)
    M = arrayio.read(args.input)
    R = jmath.start_R()
    gene_id = M._row_names[in_id]
    jmath.R_equals_vector(gene_id,'gene_id')
    R('library(biomaRt)')
    jmath.R_equals(in_attribute,'in_attribute')
    jmath.R_equals(in_attribute,'filters')
    jmath.R_equals(in_mart,'in_mart')
    R('old=useMart("ensembl",in_mart)')
    id_pair = []
    for out_platform in out_platforms:
        out_attribute = arrayplatformlib.get_bm_attribute(out_platform)
        out_mart = arrayplatformlib.get_bm_organism(out_platform)
        jmath.R_equals(out_attribute,'out_attribute')
        jmath.R_equals(out_mart,'out_mart')
        R('new=useMart("ensembl",out_mart)')
        R('homolog = getLDS(attributes=in_attribute,filters=filters,values=gene_id,mart=old,attributesL=out_attribute,martL=new)')
        homolog=R['homolog']
        old_id = [str(i) for i in homolog[0]]
        human_id = [str(i) for i in homolog[1]]
        id_pair.append((old_id,human_id))
    #align different new platform id and the matrix row name
    index_list = []
    new_id = []
    for j in range(len(id_pair)):
            new_id.append([])
    for i in range(len(M._row_names[in_id])):
        gene_name = M._row_names[in_id][i]
        num_index = []
        gene_index_list = []
        for j in range(len(id_pair)):
            gene_index = [index for index,item in enumerate(id_pair[j][0]) if item == str(gene_name)]
            num_index.append(len(gene_index))
            gene_index_list.append(gene_index)
        max_n = max(num_index)
        if max_n == 0:
            index_list.append(i)
            for j in range(len(id_pair)): 
                new_id[j].append('')
        else:
            index_list.extend([M._row_names[in_id].index(gene_name)]*max_n)
            for h in range(len(id_pair)):
                if len(gene_index_list[h])<max_n:
                    gene_index_list[h].extend([gene_index_list[h][0]]*(max_n-len(gene_index_list[h])))
                a=[id_pair[h][1][k] for k in gene_index_list[h]]
                new_id[h].extend(a)
    M = M.matrix(index_list,None)
    ids = M._row_order
    id_len = len(ids)
    #change the platform name if we will add the exact same headers
    ids.extend(out_platforms)
    for i in range(len(out_platforms)):
        name_index = [index for index,item in enumerate(ids) if item == out_platforms[i]]
        for j in range(1,len(name_index)):
            ids[name_index[j]] = ids[name_index[j]]+'_'+str(j)
            out_platforms[name_index[j]-id_len] = out_platforms[name_index[j]-id_len]+'_'+str(j)
    M._row_order=ids
    for i in range(len(out_platforms)):
        M._row_names[out_platforms[i]]=new_id[i]
    newfile=file(args.outpath,'w')
    arrayio.tab_delimited_format.write(M,newfile)
    newfile.close()
    os.chdir(cwd)

            
if __name__=='__main__':
    main()
