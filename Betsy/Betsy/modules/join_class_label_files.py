#join_class_label_file.py
import read_label_file
import module_utils
import shutil
import os

def run(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    clf1 = module_utils.find_object(
        parameters,objects,'class_label_file','merge1')
    clf2 = module_utils.find_object(
        parameters,objects,'class_label_file','merge2')
    assert os.path.exists(clf1.identifier),(
        'the input clf1 %s for join_class_label_file does not exist'
        %clf1.identifier)
    assert os.path.exists(clf2.identifier),(
        'the input clf2 %s for join_class_label_file does not exist'
        %clf2.identifier)
    result1,label_line1,second_line1=read_label_file.read(clf1.identifier)
    result2,label_line2,second_line2=read_label_file.read(clf2.identifier)
    second_line1.extend(second_line2)
    label = [str(int(i)+len(result1)) for i in label_line2]
    label_line1.extend(label)
    read_label_file.write(outfile,second_line1,label_line1)
    assert module_utils.exists_nz(outfile),(
        'the output file %s for join_class_label_file fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
                                  parameters,[clf1,clf2],pipeline,outfile)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    parameters = module_utils.renew_parameters(
             parameters,['merge1','merge2'])
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)
    
def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'class_label_file_'+original_file+'.cls'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'class_label_file','merge1')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for join_class_label_files does not exist'
        %single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    parameters = module_utils.renew_parameters(
             parameters,['merge1','merge2'])
    attributes = parameters.values()
    new_object = module_utils.DataObject('class_label_file',attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects

    
