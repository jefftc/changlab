#convert_xls_tdf.py
import os
import module_utils
import xlrd

def run(parameters,objects):
    """convert xls signal file to pcl format"""
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects)
    #convert it to tdf format
    f = file(outfile,'w')
    book = xlrd.open_workbook(identifier)
    sheet = book.sheet_by_index(0)
    for i in range(sheet.nrows):
        for j in range(sheet.ncols):
            c=sheet.cell(i,j).value
            f.write(str(c)+'\t')
        f.write('\n')
    f.close()
    #convert to pcl format and rewrite 
    import arrayio
    M = arrayio.read(outfile)
    M_c = arrayio.convert(M,to_format=arrayio.pcl_format)
    f = file(outfile,'w')
    arrayio.pcl_format.write(M_c,f)
    f.close()
    module_utils.write_Betsy_parameters_file(parameters,single_object)
    return new_objects


def make_unique_hash(parameters,objects):
    return module_utils.make_unique_hash(parameters,objects,'signal_file','Contents,DatasetId')

def get_outfile(parameters,objects):
    return module_utils.get_outfile(parameters,objects,'signal_file','Contents,DatasetId','signal_file')
    
def get_identifier(parameters,objects):
    return module_utils.find_object(parameters,objects,'signal_file','Contents,DatasetId')
