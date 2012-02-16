#module for reading class label files
#read_label_file.py

def find_index_list(inputlist,key):
    """get a list of index for key in inputlist"""
    start = 0
    indexlist = []
    while 1:
        try:
            index = inputlist.index(key,start)
        except:
             break
        indexlist.append(index)
        start = index+1
    return indexlist

def read(filename):
    """read the cls file and return a list of tuples(label_index,label)"""
    f=open(filename,'r')
    content=f.readlines()
    f.close()
    assert content[0] != '#numeric','numeric label file'
    first_line = content[0].split()
    assert first_line[2] != 1, 'file is not cls format'
    sample_num = int(first_line[0])
    class_num = int(first_line[1])
    second_line = content[1]
    assert second_line.startswith('#'),'file is not cls format'
    second_line = second_line.split()# ['#', 'ALL', 'AML']
    assert len(second_line)-1 == class_num, 'number of class is not consistent'
    label_line= content[2].strip().split()
    assert len(label_line) ==  sample_num,' number of sample is not consistent'
    label_dict = dict()
    for i in label_line:
        label_dict[i] = label_dict.get(i,0) + 1
    result = []
    for i in range(class_num):
         single_class = find_index_list(label_line,label_dict.keys()[i])
         result.append((single_class,label_dict.keys()[i]))    
    return result,label_line,second_line[1:]

def write(filename,class_name,label_line):
    f = file(filename,'w')
    f.write('\t'.join([str(len(label_line)),str(len(class_name)),'1']))
    f.write('\n')
    a=['#']
    a.extend(class_name)
    f.write('\t'.join(a))
    f.write('\n')
    f.write('\t'.join(label_line))
    f.close()
