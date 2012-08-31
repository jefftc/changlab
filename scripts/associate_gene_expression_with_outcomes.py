#!/usr/bin/env python

from genomicode import filelib,jmath,genesetlib,config,matrixlib
import arrayio
import os
import argparse
import numpy

def main():

    parser = argparse.ArgumentParser(
        description = 'associate the gene expression with outcomes')
    parser.add_argument('-o',dest = 'outpath',type = str,
                        help = 'specify the outpath',default=None)
    parser.add_argument('-f',dest = 'expression_file',type = str,
                        help = 'specify the gene expression file',default=None)
    parser.add_argument('-cf',dest = 'clinical_data',type = str,
                        help = 'specify the clinical data file', default=None)
    parser.add_argument('--time_header',dest = 'time_header',type = str,
                        help = 'specify the header for time column', default=None)
    parser.add_argument('--dead_header',dest = 'dead_header',type = str,
                        help = 'specify header for dead column', default=None)
    parser.add_argument('--genes',dest = 'genes',type = str,
                        help = 'specify genes for analysis', default=None,
                        action='append')
    parser.add_argument('--prism_file',dest = 'prism_file',type = str,
                        help = 'specify the prism file name', default=None)
    args = parser.parse_args()
    outfile = args.outpath
    assert outfile,'please specify the path of output file'
    input_file = args.expression_file
    assert input_file,'please specify the path of gene expression data file'
    clinical_file = args.clinical_data
    assert clinical_file,'please specify the path of clinical data'
    time_header = args.time_header
    assert time_header,'please specify the header of month column in clinical data'
    dead_header = args.dead_header
    assert dead_header,'please specify the header of dead column in clinical data'
    genes = args.genes
    if genes:
        gene_list=[]
        for i in genes:
            gene_list.extend(i.split(','))
            
    clinical_data = genesetlib.read_tdf(clinical_file,preserve_spaces=True, allow_duplicates=True)
    time_data = None
    dead_data = None
    name_in_clinical = None
    M = arrayio.read(input_file)
    sample_name = M._col_names['_SAMPLE_NAME']
       
    for g in clinical_data:
        n = len(list(set(sample_name).intersection(set(g[2]))))
        if n > 0:
           name_in_clinical = g[2]
        if g[0] == time_header:
            time_data = g[2]
        if g[0] == dead_header:
            dead_data = g[2]
    assert name_in_clinical,'cannot match the sample name in clinical data and gene expression data'
    assert time_data,'there is no column named %s'%time_header
    assert dead_data,'there is no column named %s'%dead_header
    #only consider the sample who has month and dead information
    sample_index1 = [index for index,item in enumerate(time_data) if len(item)>0]
    sample_index2 = [index for index,item in enumerate(dead_data) if len(item)>0]
    sample_index = list(set(sample_index1).intersection(set(sample_index2)))
    
    #align the sample order of gene expression file and clinical data
    x= matrixlib.align_cols_to_annot(M,name_in_clinical,reorder_MATRIX=True)
    M,colnames = x
    #if given genes,get the row index of match genes
    if genes:
        row_index = []
        for gene in gene_list:
            for column_name in M._row_order:
                if gene in M.row_names(column_name):
                    index_list = [index for index,item in enumerate(M.row_names(column_name)) if item == gene]
                    row_index.extend(index_list)
        assert row_index,'we cannot find the genes you given in the expression data'
        M = M.matrix(row_index,None)
        
    ids = M._row_order
    geneid = M._row_names[ids[0]]
    data_all = M.slice()
    output = []
    R=jmath.start_R()
    for i in range(len(geneid)):
        data = data_all[i]
        data_new = [data[j] for j in sample_index]
        median = numpy.median(data_new)
        group1= [index for index,item in enumerate(data) if item>=median and index in sample_index]
        group2= [index for index,item in enumerate(data) if item<median and index in sample_index]
        group1_month=[float(time_data[k]) for k in group1]
        group2_month=[float(time_data[k]) for k in group2]
        survival1=numpy.mean(group1)
        survival2=numpy.mean(group2)
        group1_dead=[int(dead_data[k]) for k in group1]
        group2_dead=[int(dead_data[k]) for k in group2]
        med_high = numpy.median(group1_month)
        med_low = numpy.median(group2_month)
        jmath.R_equals(group1_month,'group1')
        jmath.R_equals(group2_month,'group2')
        jmath.R_equals(group1_dead,'dead1')
        jmath.R_equals(group2_dead,'dead2')
        R('source("'+config.kaplanmeierlib+'")')
        R('a<-calc.km(group1,dead1,group2,dead2)')
        c=R['a']
        if c[0][0]<0.05:
            if survival1>survival2:
                direction = 'High gene expression correlates with high survival.'
            else:
                direction = 'High gene expression correlates with low survival.'
        else:
            direction = ''
        output.append([str(c[0][0]),str(len(group2_month)),
                       str(len(group1_month)),str(med_low),str(med_high),direction,str(median)])
        if args.prism_file:
            assert len(geneid)==1,'multiple genes match and cannot write the prism file'
            jmath.R_equals('"'+args.prism_file+'"','filename')
            R('write.km.prism(filename,group1,dead1,group2,dead2)') 
    f=file(outfile,'w')
    headers=ids[:]
    headers.extend(['p.50%','# patients low','# patients high','median survival low',
                    'median survival high','Direction','Cutoff'])
    f.write('\t'.join(headers)+'\n')
    for i in range(len(geneid)):
       descriptions = [M._row_names[name][i] for name in ids]
       descriptions.extend(output[i])
       f.write('\t'.join(descriptions)+'\n')
    f.close()
               
if __name__=='__main__':
    main()
