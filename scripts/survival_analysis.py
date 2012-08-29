#!/usr/bin/env python

from genomicode import filelib,jmath,genesetlib,config
import arrayio
import os
import argparse
import numpy
def gunzip(filename):
    import gzip
    if filename.endswith('.gz'):
        newfilename =os.path.join(os.getcwd(),os.path.split(os.path.splitext(filename)[0])[-1])
        #unzip the gz data
        fileObj = gzip.GzipFile(filename, 'rb');
        fileObjOut = file(newfilename, 'wb');
        while 1:
            line = fileObj.readline()
            if line == '':
                break
            fileObjOut.write(line)
        fileObj.close()
        fileObjOut.close()
        assert os.path.exists(newfilename),'unzip the file %s fails'%filename
        return newfilename
    
def _transpose_gmx(matrix):
    # GMX format:
    # <gene set name>  ...
    # <description>    ...
    # <gene>           ...
    #
    # Each column is a gene set.
    t_matrix = []
    # Iterate over each of the columns (gene sets).
    for j in range(len(matrix[0])):
        x = []
        # Pull this column out of the matrix into variable x.
        for i in range(len(matrix)):
            # These lines may not be the same length as the names,
            # e.g. if the gene sets at the end have very few genes.
            if j >= len(matrix[i]):
                continue
            x.append(matrix[i][j])
        assert len(x) >= 2
        x1 = x[:2]
        x2 = [x.strip() for x in x[2:]]
        x = x1 + x2
        t_matrix.append(x)
    # The rows are not guaranteed to be the same length.
    return t_matrix
def read_tdf(filename):
    # yield name, description (always ""), list of genes
    #import filelib
    matrix = [x for x in filelib.read_cols(filename)]
    t_matrix = _transpose_gmx(matrix)
    # For the TDF file, each of the rows should be exactly the same
    # length.  Make sure they are the same length.
    maxlen = max([len(x) for x in t_matrix]) - 1
    assert maxlen >= 0
    cli_dict = {}
    for i in range(len(t_matrix)):
        x = t_matrix[i]
        name, genes = x[0],  x[1:]
        cli_dict[name] = genes
    return cli_dict

def main():

    parser = argparse.ArgumentParser(
        description = 'survival analysis')
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
    args = parser.parse_args()
    cwd = os.getcwd()
    outfile = args.outpath
    assert outfile,'please specify the path of output file'
    input_file = args.expression_file
    if input_file.endswith('gz'):
        input_file = gunzip(input_file)
    assert input_file,'please specify the path of gene expression data file'
    clinical_file = args.clinical_data
    assert clinical_file,'please specify the path of clinical data'
    time_header = args.time_header
    assert time_header,'please specify the header of month column in clinical data'
    dead_header = args.dead_header
    assert dead_header,'please specify the header of dead column in clinical data'
    
    clinical_dict = read_tdf(clinical_file)
    #only consider the sample who has month and dead information
    sample_index = [index for index,item in enumerate(clinical_dict[time_header]) if len(item)>0]

    M = arrayio.read(input_file)
    #check the sample order of gene expression file and clinical data
    assert M._col_names['_SAMPLE_NAME'] in clinical_dict.values(),(
        'the order of sample in gene_expression file and clinical data is different')
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
        group1_month=[float(clinical_dict[time_header][i]) for i in group1]
        group2_month=[float(clinical_dict[time_header][i]) for i in group2]
        group1_dead=[int(clinical_dict[dead_header][i]) for i in group1]
        group2_dead=[int(clinical_dict[dead_header][i]) for i in group2]
        med_high = numpy.median(group1_month)
        med_low = numpy.median(group2_month)
        jmath.R_equals(group1_month,'group1')
        jmath.R_equals(group2_month,'group2')
        jmath.R_equals(group1_dead,'dead1')
        jmath.R_equals(group2_dead,'dead2')
        R('source("'+config.kaplanmeierlib+'")')
        R('a<-calc.km(group1,dead1,group2,dead2)')
        c=R['a']
        output.append([str(c[0][0]),str(len(group2_month)),str(len(group1_month)),str(med_low),str(med_high)])

    f=file(outfile,'w')
    headers=ids[:]
    headers.extend(['p.50%','# patients low','# patients high','median survival low', 'median survival high'])
    f.write('\t'.join(headers)+'\n')
    for i in range(len(geneid)):
       descriptions = [M._row_names[name][i] for name in ids]
       descriptions.extend(output[i])
       f.write('\t'.join(descriptions)+'\n')
    f.close()
    os.chdir(cwd)

            
if __name__=='__main__':
    main()
