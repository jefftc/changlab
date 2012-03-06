#plot_hyb_bar.py
import os
import module_utils
import shutil

def run(parameters,objects,pipeline):
    identifier,single_object = get_identifier(parameters,objects)
    outfile,new_objects = get_outfile(parameters,objects,pipeline)
    plot_hyb_bar(identifier,outfile)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline)
    return new_objects

def make_unique_hash(parameters,objects,pipeline):
    return module_utils.make_unique_hash(
        parameters,objects,'signal_file','Contents,DatasetId',pipeline)

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'signal_file','Contents,DatasetId','plot_hyb_bar',pipeline)
    
def get_identifier(parameters,objects):
    return module_utils.find_object(
        parameters,objects,'signal_file','Contents,DatasetId')

def plot_hyb_bar(filename,outfile):
    high = ['ILMN_2038770','ILMN_2038769']
    med = ['ILMN_2038768','ILMN_2038771']
    low = ['ILMN_1343050','ILMN_1343052']
    high_data = []
    med_data = []
    low_data = []
    from genomicode import jmath
    import arrayio
    M = arrayio.read(filename)
    header = M.row_names()
    for i in range(M.dim()[0]):
        if not M.row_names(header[1])[i] == 'cy3_hyb':
            continue
        if M.row_names(header[0])[i] in high:
            high_data.extend(M.slice()[i])
        if M.row_names(header[0])[i] in med:
            med_data.extend(M.slice()[i])
        if M.row_names(header[0])[i] in low:
            low_data.extend(M.slice()[i])
    high_data.extend(med_data)
    high_data.extend(low_data)
    row = len(high_data)/3.0
    assert row > 0,'input is not a control file'
    R = jmath.start_R()
    jmath.R_equals_vector(high_data,'y')
    jmath.R_equals(row,'row')
    R(' y <- matrix(y,row,3)')
    R(' y.means <- apply(y,2,mean)')
    R('y.sd <- apply(y,2,sd)')
    R('library(R.utils)')
    command = 'png2("' + outfile + '")'
    R(command)
    R('barx <- barplot(y.means, names.arg=c("high","med","low"), \
      col="blue", axis.lty=1, ylab="Signal")')
    R('arrows(barx,y.means+1.96*y.sd/10, barx, y.means-1.96*y.sd/10, \
      angle=90, code=3, length=0.1)')
    R('dev.off()')
    assert module_utils.exists_nz(outfile),'the plot_hyb_bar.py fails'
        
