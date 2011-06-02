"""This program is to apply a log2 algorithm on the gct file
    Output file is a gct file with 'L2' append in the filename"""
import os
import sys
import math

def main(filename):
    
    a = os.getcwd()
    assert os.path.exists(os.path.join(a,filename)),'the file does not exist'
    handle = open(filename)
    line = handle.read().split('\n')
    
    outfilename = filename[:-4]+'L2'+filename[-4:]
    f = open(outfilename,'w')
    
    #write the first three line of gct file
    for i in range(3):
        f.write(line[i]+'\n')
        
    #log the data and write it to the file
    for i in range(3,len(line)):
        number = line[i].split('\t')
        if len(number)>1:
            newnumber = [number[0],number[1]]
            for k in range(2,len(number)):
                a = math.log(float(number[k]),2)
                newnumber.append(str(a))
            f.write("\t".join(newnumber)+"\n")
     
    print 'yes,log'
    f.close()


if __name__=='__main__':

    main(sys.argv[1])
