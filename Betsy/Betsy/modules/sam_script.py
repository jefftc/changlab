#sam_script.py

import shutil
import os
import arrayio
from Betsy import read_label_file
from genomicode import jmath
import sys

def main(filename,label_file,outfile,delta,foldchange):
    label,label_line,second_line = read_label_file.read(label_file)
    M = arrayio.read(filename)
    data = M.slice()
    label_list = [int(i)+1 for i in label_line]
    key = M._row_names.keys()
    genenames = M._row_names[key[0]]
    genenames = [str('\"'+i+'\"') for i in genenames]
    if not os.path.exists(outfile):
        os.mkdir(outfile)
    pngfig = os.path.join(outfile,'sam_plot.png')
    out_result = os.path.join(outfile,'sam_reult.txt')
    gene_ids,result = sam(data,label_list,genenames,delta,foldchange,pngfig)
    group1s = []
    group2s = []
    sd1s = []
    sd2s = []
    label1 = second_line[int(label[0][1])]
    label2 = second_line[int(label[1][1])]
    if gene_ids:
        index_list = []
        for key in M._row_order:
            full_gene_list = M._row_names[key]
            if gene_ids[0] not in full_gene_list:
                continue
            index_list = [full_gene_list.index(i) for i in gene_ids]
        M_select = M.slice(index_list,None)
        first = []
        second = []
        for i in range(len(M_select)):
             first.append([M_select[i][j] for j in label[0][0]])
             second.append([M_select[i][j] for j in label[1][0]])
        higher_group=[]
        for i in range(len(first)):
                group1 = sum(first[i])/float(len(first[i]))
                group2 = sum(second[i])/float(len(second[i]))
                sd1 = jmath.stddev_list(first[i])
                sd2 = jmath.stddev_list(second[i])
                group1s.append(group1)
                group2s.append(group2)
                sd1s.append(sd1)
                sd2s.append(sd2)
                if group1 >= group2:
                    higher_group.append(label1)
                else:
                    higher_group.append(label2)
                    
    header = ['Gene Name','Ave_'+label1,
          'Ave_'+label2,'SD_'+label1,'SD_'+label2,'Score(d)',
          'Numerator(r)','Denominator(s+s0)',
          'Fold Change','q_value','higher_expression']
    f=file(out_result,'w')
    f.write('\t'.join(header))
    f.write('\n')
    if gene_ids:
        for i in range(len(gene_ids)):
            f.write(str(gene_ids[i])+'\t')
            f.write(str(group1s[i])+'\t')
            f.write(str(group2s[i])+'\t')
            f.write(str(sd1s[i])+'\t')
            f.write(str(sd2s[i])+'\t')
            for k in range(5):
                f.write(str(result[k][i])+'\t')
            f.write(str(higher_group[i])+'\n')
    f.close()
    

   
def sam(X, Y, genenames,delta,foldchange,pngfig):
    """X is a matrix slice,Y is the label list"""
    assert len(X[0]) == len(Y), 'X and Y should be equal length'
    R = jmath.start_R()
    jmath.R_equals_matrix(X,'x',True)
    jmath.R_equals(genenames,'genenames')
    jmath.R_equals(Y,'y')
    jmath.R_equals(foldchange,'foldchange')
    jmath.R_equals(delta,'DELTA')
    R('library(samr)')
    R('D<-list(x=x,y=y,logged2=TRUE,geneid = 1:length(x),genenames=genenames)')
    R('S<-samr(D,resp.type="Two class unpaired",nperms=100)')
    R('DTAB<-samr.compute.delta.table(S,min.foldchange=foldchange)')
    R('SIG<-samr.compute.siggenes.table(S,DELTA,D,DTAB,min.foldchange=foldchange)')
    R('up<-SIG$ngenes.up')
    R('lo<-SIG$ngenes.lo')
    R('library(R.utils)')
    command = 'png2("'+pngfig+'")'
    R(command)
    R('samr.plot(S,DELTA)')
    R('title(main=paste("DELTA=",DELTA))')
    R('dev.off()')
    import rpy2.robjects as robjects
    R = robjects.r
    up = R['up']
    lo = R['lo']
    gene_ids = []
    scores = []
    numerators = []
    denominators = []
    foldchanges = []
    q_values = []
    if up[0] > 0:
        R('geneID1<-SIG$genes.up[,"Gene ID"]')
        R('Score1<-SIG$genes.up[,"Score(d)"]')
        R('Numerator1<-SIG$genes.up[,"Numerator(r)"]')
        R('Denominator1<-SIG$genes.up[,"Denominator(s+s0)"]')
        R('foldchange1<-SIG$genes.up[,"Fold Change"]')
        R('q1<-SIG$genes.up[,"q-value(%)"]')
        gene_id=R['geneID1']
        gene_ids.extend(gene_id)
        scores.extend(R['Score1'])
        numerators.extend(R['Numerator1'])
        denominators.extend(R['Denominator1'])
        foldchanges.extend(R['foldchange1'])
        q = [float(i)/100 for i in R['q1']]
        q_values.extend(q)
        
    if lo[0] > 0:
        R('geneID2<-SIG$genes.lo[,"Gene ID"]')
        R('Score2<-SIG$genes.lo[,"Score(d)"]')
        R('Numerator2<-SIG$genes.lo[,"Numerator(r)"]')
        R('Denominator2<-SIG$genes.lo[,"Denominator(s+s0)"]')
        R('foldchange2<-SIG$genes.lo[,"Fold Change"]')
        R('q2<-SIG$genes.lo[,"q-value(%)"]')
        gene_id=R['geneID2']
        gene_ids.extend(gene_id)
        scores.extend(R['Score2'])
        numerators.extend(R['Numerator2'])
        denominators.extend(R['Denominator2'])
        foldchanges.extend(R['foldchange2'])
        q = [float(i)/100 for i in R['q2']]
        q_values.extend(q)
        
    return gene_ids,[scores,numerators,denominators,foldchanges,q_values]

if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2],sys.argv[3],float(sys.argv[4]),float(sys.argv[5]))
