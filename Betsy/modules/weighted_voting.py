#weighted_voting.py
import module_utils
import shutil
import os
from genomicode import jmath
import Betsy_config
import subprocess
import read_label_file
def run(parameters,objects,pipeline):
    train_identifier,single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    train_label_file,obj = module_utils.find_object(
        parameters,objects,'class_label_file','TrainContents,DatasetId')
    test_label_file,obj = module_utils.find_object(
        parameters,objects,'class_label_file','TestContents,DatasetId')
    test_file,obj = module_utils.find_object(
        parameters,objects,'signal_file','TestContents,DatasetId')
    assert os.path.exists(test_file),'the test file does not exist'
    assert os.path.exists(test_label_file),'cannot find test_label_file'
    assert os.path.exists(train_label_file),'cannot find train_label_file'
    module_name = 'WeightedVoting'
    gp_parameters = dict()
    gp_parameters['train.filename'] = train_identifier
    gp_parameters['train.class.filename'] = train_label_file
    gp_parameters['test.filename'] = test_file
    gp_parameters['test.class.filename'] = test_label_file
    gp_module = Betsy_config.GENEPATTERN
    command = [gp_module, module_name]
    for key in gp_parameters.keys():
        a = ['--parameters',key+':'+ gp_parameters[key]]
        command.extend(a)
    
    download_directory = None
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    process.wait()
    output_text =  process.stdout.read()
    out_lines = output_text.split('\n')
    for out_line in out_lines:
        if out_line != 'Loading required package: rJava' and len(out_line)>0:
            download_directory = out_line
            break
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    result_files = os.listdir(download_directory)
    assert 'stderr.txt' not in result_files,'gene_pattern get error'
    gp_files = os.listdir(download_directory)
    result,label_line,second_line = read_label_file.read(train_label_file)
    for gp_file in gp_files:
        if gp_file.endswith('pred.odf'):
            gp_file = os.path.join(download_directory,gp_file)
            f = file(gp_file,'r')
            text = f.readlines()
            assert text[1][0:12]=='HeaderLines='
            start=int(text[1][12:-1])
            label=[]
            for i in text[start+2:]:
                line = i.split()
                label.append(line[2])
            number = [second_line.index(i) for i in label]
            legend_name = ['"'+ i +'='+
                           str(second_line.index(i))+'"' for i in second_line]
            R = jmath.start_R()
            R('library(R.utils)')
            command = 'png2("'+outfile+'")'
            R(command)
            jmath.R_equals_vector(number,'p_label')
            jmath.R_equals_vector(legend_name,'legend_name')
            R('barx <- barplot(p_label, ylim=c(-1,max(p_label)+1), \
              col="blue", axis.lty=2, ylab="Prediction",xlab="Sample")')
            R('legend("bottomleft", legend_name, lty=1,cex=1)')
            R('dev.off()')
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'signal_file','TrainContents,DatasetId')
    assert os.path.exists(identifier),'the train file does not exist'
    return identifier,single_object

def get_outfile(parameters,objects,pipeline):
    return module_utils.get_outfile(
        parameters,objects,'signal_file','TrainContents,DatasetId',pipeline)
    
def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'weighted_voting',parameters,objects,single_object)
    return new_objects
