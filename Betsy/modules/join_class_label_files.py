#join_class_label_file.py
import read_label_file
import module_utils
import shutil
import rule_engine
import os

def run(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    clf1,obj1 = module_utils.find_object(
        parameters,objects,'class_label_file','merge1,dataset1')
    clf2,obj2 = module_utils.find_object(
        parameters,objects,'class_label_file','merge2,dataset2')
    assert os.path.exists(clf1),'the input clf1 does not exist'
    assert os.path.exists(clf2),'the input clf2 does not exist'
    result1,label_line1,second_line1=read_label_file.read(clf1)
    result2,label_line2,second_line2=read_label_file.read(clf2)
    second_line1.extend(second_line2)
    label = [str(int(i)+len(result1)) for i in label_line2]
    label_line1.extend(label)
    read_label_file.write(outfile,second_line1,label_line1)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(parameters,[obj1,obj2],pipeline)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    parameters = module_utils.renew_parameters(
             parameters,['merge1','merge2','dataset1','dataset2'])
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)
    
def get_outfile(parameters,objects,pipeline):
    outfile = os.path.join(os.getcwd(),'class_label_file.cls')
    return outfile

def get_identifier(parameters,objects):
    identifier,single_object = module_utils.find_object(
        parameters,objects,'class_label_file','merge1,dataset1')
    return identifier,single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    identifier,single_object = get_identifier(parameters,objects)
    parameters = module_utils.renew_parameters(
             parameters,['merge1','merge2','dataset1','dataset2'])
    attributes = parameters.values()
    new_object = rule_engine.DataObject('class_label_file',attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects

    
