#!/usr/bin/env python

from genomicode import filelib,jmath,genesetlib,config,matrixlib
import arrayio
import os
import argparse
import numpy
import sys
def main():
    parser = argparse.ArgumentParser(
        description = 'associate the gene expression with outcomes')
    parser.add_argument(dest = 'expression_file',type = str,
                        help = 'specify the gene expression data file', default=None)
    parser.add_argument(dest = 'clinical_data',type = str,
                        help = 'specify the clinical data file', default=None)
    parser.add_argument('--outcome',dest = 'outcome',type = str,
                        help = 'specify time_header and dead_header in the format<time_header>,<dead_header>',default=[],
                        action='append')
    parser.add_argument('--genes',dest = 'genes',type = str,
                        help = 'specify genes for analysis', default=None,
                        action='append')
    parser.add_argument('-o',dest = 'outpath',type = str,
                        help = 'specify the outpath',default=None)
    parser.add_argument('--prism_file',dest = 'prism_file',type = str,
                        help = 'specify the prism file name', default=None)
    parser.add_argument('--cutoff',dest = 'cutoff',type = float,
                        help = 'specify the cutoff(between 0 and 1) to calculate',
                        default=[],action='append')
    args = parser.parse_args()

    input_file = args.expression_file
    assert input_file,'please specify the path of gene expression data file'
    clinical_file = args.clinical_data
    assert clinical_file,'please specify the path of clinical data'
    outcomes = args.outcome
    assert len(outcomes)>0,'please specify the time_header and dead_header'
    cutoffs = args.cutoff
    if not cutoffs:
        cutoffs = [0.5]
    for cutoff in cutoffs:
        assert cutoff>0 and cutoff<1.0, 'cutoff should be between 0 and 1'
    genes = args.genes
    if genes:
        gene_list=[]
        for i in genes:
            gene_list.extend(i.split(','))
            
    M = arrayio.read(input_file)
    clinical_data = genesetlib.read_tdf(clinical_file,preserve_spaces=True, allow_duplicates=True)
    sample_name = M._col_names['_SAMPLE_NAME']
    name_in_clinical = None
    clinical_dict = dict()
    for g in clinical_data:
        clinical_dict[g[0]] = g[2]
    for g in clinical_dict:
        n = len(list(set(sample_name).intersection(set(clinical_dict[g]))))
        if n > 0:
           name_in_clinical = clinical_dict[g]
    assert name_in_clinical,'cannot match the sample name in clinical data and gene expression data'
    
    #align the sample order of gene expression file and clinical data
    x = matrixlib.align_cols_to_annot(M,name_in_clinical,reorder_MATRIX=True)
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
    headers=ids[:]
    data_all = M.slice()
    output_data = []
    for i in range(len(geneid)):
        output_data.append([M._row_names[name][i] for name in ids])
    kaplanmeierlib_path = config.kaplanmeierlib
    assert os.path.exists(kaplanmeierlib_path),(
        'can not find the kaplanmeierlib script %s'%kaplanmeierlib_path)
    R = jmath.start_R()
    R('require(splines,quietly=TRUE)')
    R('source("'+config.kaplanmeierlib+'")')
    R('ow<-options("warn")')
    for outcome in outcomes:
        time_data = None
        dead_data = None
        x = outcome.split(',')
        assert len(x)==2, 'the outcome format should be <time_header>,<dead_header>'
        time_header,dead_header = x
        for g in clinical_dict:
            if g == time_header:
                time_data = clinical_dict[g]
            if g == dead_header:
                dead_data = clinical_dict[g]
        assert time_data,'there is no column named %s'%time_header
        assert dead_data,'there is no column named %s'%dead_header
        #only consider the sample who has month and dead information
        sample_index1 = [index for index,item in enumerate(time_data) if len(item)>0]
        sample_index2 = [index for index,item in enumerate(dead_data) if len(item)>0]
        sample_index = list(set(sample_index1).intersection(set(sample_index2)))
        for percentage in cutoffs:
            outer = '(Outer '+ str(int(percentage*100)) +'%)'
            newheaders = ['p '+outer,'Num Patients Low '+outer,'Num Patients High '+outer,
                        '50% Survival Low '+outer,'50% Survival High '+outer,
                        '90% Survival Low '+outer,'90% Survival High '+outer,
                        'Relation '+outer,'Low Expression '+outer,'High Expression '+outer]
            if len(outcomes)>1:
                    newheaders = [time_header+' ' + i for i in newheaders]
            headers.extend(newheaders)
        for i in range(len(geneid)):
            data = data_all[i]
            data_new = [data[j] for j in sample_index]
            data_order = data_new[:]
            data_order.sort()
            for percentage in cutoffs:
                high_point = data_order[int(round((len(data_order)-1)*(1-percentage)))]
                low_point = data_order[int(round((len(data_order)-1)*percentage))]
                group1_index = [index for index,item in enumerate(data) if item>=high_point and index in sample_index]
                group2_index = [index for index,item in enumerate(data) if item<low_point and index in sample_index]
                assert len(group1_index)>0,'there is no patient in the high expression group for cutoff%f'%percentage
                assert len(group2_index)>0,'there is no patient in the low expression group for cutoff%f'%percentage
                group1_month=[float(time_data[k]) for k in group1_index]
                group2_month=[float(time_data[k]) for k in group2_index]
                group1_dead=[int(dead_data[k]) for k in group1_index]
                group2_dead=[int(dead_data[k]) for k in group2_index]
                jmath.R_equals(group1_month,'group1')
                jmath.R_equals(group2_month,'group2')
                jmath.R_equals(group1_dead,'dead1')
                jmath.R_equals(group2_dead,'dead2')
                R('options(warn=-1)')
                R('a<-calc.km(group1,dead1,group2,dead2)')
                R('options(ow)')
                c=R['a']
                med_high_50 = c.rx2('surv1.50')[0]
                med_high_90 = c.rx2('surv1.90')[0]
                med_low_50 = c.rx2('surv2.50')[0]
                med_low_90 = c.rx2('surv2.90')[0]
                p_value = c.rx2('p.value')[0]
                if p_value<0.05:
                    if 'NA' not in [str(med_high_50),str(med_low_50)]:
                        if med_high_50>med_low_50:
                            direction = 'High gene expression correlates with high survival.'
                        else:
                            direction = 'High gene expression correlates with low survival.'
                    elif  'NA' not in [str(med_high_90),str(med_low_90)]:
                        
                        if med_high_90>med_low_90:
                            direction = 'High gene expression correlates with high survival.'
                        else:
                            direction = 'High gene expression correlates with low survival.'
                    else:
                        direction = ''    
                else:
                    direction = ''
                output_data[i].extend([str(p_value),str(len(group2_month)),
                        str(len(group1_month)),str(med_low_50),str(med_high_50),str(med_low_90),
                        str(med_high_90),direction,str(low_point),str(high_point)])
                
            if args.prism_file:
                assert len(geneid)==1,'multiple genes match and cannot write the prism file'
                jmath.R_equals('"'+args.prism_file+'"','filename')
                R('write.km.prism(filename,"High Expression",group1,dead1,"Low Expression",group2,dead2)')
    
    f = sys.stdout
    if args.outpath: 
        f = file(args.outpath,'w')
    print >>f,'\t'.join(headers) 
    for i in range(len(geneid)):
       print >>f,'\t'.join(output_data[i])
    f.close()
               
if __name__=='__main__':
    main()
