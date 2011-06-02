"""this program is to convert a normalized file in text format to gct format"""

import os
import sys

def main(filename):
    a=os.getcwd()
    assert os.path.exists(os.path.join(a,filename)),'the file does not exist'
    handle = open(filename)
       
    line = handle.readline()
    columns = line.rstrip("\n").split("\t")
    sample_names = columns[1:len(columns)]
       
    probe_set_ids = []
    matrix = []
    while True:
           line = handle.readline().rstrip()
           if not line:
               break
           columns = line.split("\t")
       
           id = columns[0]
           values = columns[1:len(columns)]
           probe_set_ids.append(id)
           matrix.append(values)
       
    # Write out as GCT format.
    outputfile = filename[:-4]+'.gct'
    f = file(outputfile,'w')
    f.write("#1.2"+'\n')

    num_rows = len(matrix)
    num_samples = len(sample_names)
    f.write("%d\t%d" % (num_rows, num_samples)+'\n')
    x = ["NAME", "Description"] + sample_names
    f.write("\t".join(x)+'\n')
       
    i = 0
    while i < len(matrix):
           x = [probe_set_ids[i], "DESC"] + matrix[i]
           f.write("\t".join(x)+"\n")
           i = i + 1
    f.close()

if __name__=='__main__':
    main(sys.argv[1])
