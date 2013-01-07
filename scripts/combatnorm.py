#!/usr/bin/env python

import arrayio
import os
import argparse
from Betsy import read_label_file, module_utils
from genomicode import jmath, config

def main():
    
    parser = argparse.ArgumentParser(
        description = 'run the combat normalization')
    parser.add_argument('-o',dest = 'outpath',type = str,
                        help = 'specify the outpath',default='Adjusted_EIF.dat_.pcl')
    parser.add_argument('-f',dest = 'input',type = str,
                        help = 'specify the input file',default=None)
    parser.add_argument('-label',dest = 'label_file',type = str,
                        help = 'specify the input file',default=None)
    args = parser.parse_args()
    
    cwd = os.getcwd()
    filename = os.path.split(args.input)[-1]
    dirname = os.path.join(cwd,(filename +'_combat'))
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    M = arrayio.read(args.input)
    os.chdir(dirname)
    data = M.slice()
    EIF = 'EIF.dat'
    SIF = 'SIF.dat'
    #write EIF file
    header = M.col_names()[0]
    sample_names = M.col_names(header)
    f = file(EIF,'w')
    f.write('\t'.join(sample_names)+'\n')
    for line in data:
        for i in range(len(line)):
            if line[i]<1:
                line[i] = 1
        line = [str(i) for i in line]
        f.write('\t'.join(line)+'\n')
    f.close()
    #write SIF file
    a,b,c = read_label_file.read(args.label_file)
    batch = [str(int(i)+1) for i in b]
    arraynames = ['Array%d'%(i+1) for i in range(M.ncol())]
    sif_header=['Array name','Sample name','Batch']
    f = file(SIF,'w')
    f.write('\t'.join(sif_header)+'\n')
    for i in range(M.ncol()):
        f.write('\t'.join([sample_names[i],arraynames[i],batch[i]])+'\n')
    f.close()
    #run combat
    import R
    run_combat = config.run_combat
    assert os.path.exists(run_combat),'cannot find the %s' %run_combat
    R.run_R(run_combat)
    assert module_utils.exists_nz('Adjusted_EIF.dat_.xls'),('the '
                       'adjusted_EIF.dat_.xls does not exist')
    f=file('Adjusted_EIF.dat_.xls','r')
    text = f.read().split('\n')
    f.close()
    text = [i for i in text if len(i)>0]
    assert len(text) == M.nrow()+1,'the number of lines in combat result is not right'
    assert len(text[0].split('\t')) == M.ncol()+1, 'the number of column in combat result is not right'
    for i in range(M.nrow()):
        M._X[i]=text[i+1].split('\t')
    newfile=file(args.outpath,'w')
    arrayio.pcl_format.write(M,newfile)
    newfile.close()
    os.chdir(cwd)

 
if __name__=='__main__':
    main()
