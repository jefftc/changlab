#gene_ranking.py
import sys
sys.path.append('/Users/xchen2/Documents/prolog files')
import arrayio
import read_label_file

def t_test(X, Y, exact=True):
    """X,Y is a matrix slice"""
    import rpy2.robjects as robjects
    R = robjects.r
    t_value = []
    p_value = []
    assert len(X) == len(Y), 'X and Y should be equal length'
    for i in range(len(X)):
        x = robjects.FloatVector(X[i])
        y = robjects.FloatVector(Y[i])
        result = R['t.test'](x, y, exact=exact)
        t_value.append(list(result.rx2('statistic').rx2('t'))[0])
        p_value.append(list(result.rx2("p.value"))[0])
    return t_value,p_value

def correlation(M, Y):
    """M is a data matrix,Y is list of label"""
    assert M.dim()[1] == len(Y), 'the number of sample is different from labels'
    import rpy2.robjects as robjects
    R = robjects.r
    cor_value = []
    y = robjects.FloatVector(Y)
    for i in range(0,M.dim()[0]):
        x = robjects.FloatVector(M.slice()[i])
        result = R['cor'](x, y)
        cor_value.append(abs(result[0]))
    return cor_value

def find_sorted_index(inputlist,gene_list):
    """get a list of index for mapping aftersort to inputlist"""
    indexlist = []
    for key in gene_list:
        index = inputlist.index(key)
        try:
            inputlist.index(key,index+1)
        except ValueError:  
            indexlist.append(index)
    return indexlist

def correlation_for_file(data_file,label_file,gene_num=True):
    """given data_file,label_file and the number of selected
       gene,return a list of select gene name"""
    #obtain the class label
    label,label_line = read_label_file.read(label_file)
    #read the data_file and caculate the correlation
    M = arrayio.read(data_file)
    p = correlation(M,label_line)
    #sort the correlation value and obtain the list of the gene index after sorting 
    c = sorted(p,reverse=True)
    sortlist = find_sorted_index(p,c)
    #obtain the gene name in the data_file
    f = open(data_file)
    a = f.read().split('\n')
    index = 0 #for pcl file
    startrows = 2 #for pcl file
    genelist = []
    for i in range(startrows,len(a)):
        genelist.append(a[i].split('\t')[index])
    #get a list of selected gene name 
    if gene_num is not True:
        select_genelist = [genelist[i]for i in sortlist[0:gene_num]]
    else:
        select_genelist = [genelist[i]for i in sortlist]
    return select_genelist
       
def t_test_for_file(data_file,label_file,gene_num=True):
    """given data_file,label_file and the number of selected
       gene,return a list of select gene name"""
    #obtain the class label
    label,label_line,second_line = read_label_file.read(label_file)
    class_num = len(label)
    assert class_num == 2, 'the number of class is not 2'
    #read the data_file and caculate the t-test
    M = arrayio.read(data_file)
    first = M.slice(None,label[0][0])
    second = M.slice(None,label[1][0])
    t,p=t_test(first,second)
    #sort the p value and obtain the list of the gene index after sorting 
    c = sorted(p)
    sortlist = find_sorted_index(p,c)
    #obtain the gene name in the data_file
    f = open(data_file)
    a = f.read().split('\n')
    index = 0 #for pcl file
    startrows = 2 #for pcl file
    genelist = []
    for i in range(startrows,len(a)):
        genelist.append(a[i].split('\t')[index])
    #get a list of selected gene name
    if gene_num is not True:
        select_genelist = [genelist[i]for i in sortlist[0:gene_num]]
    else:
        select_genelist = [genelist[i]for i in sortlist]
    return select_genelist
    

