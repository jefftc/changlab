#join_class_label_file.py
import read_label_file
import module_utils
import shutil
import rule_engine
import os

def run(parameters,objects,pipeline):
    outfile,new_objects = get_outfile(parameters,objects)
    clf1,obj1=module_utils.find_object(parameters,objects,'class_label_file','merge1,dataset1')
    clf2,obj2=module_utils.find_object(parameters,objects,'class_label_file','merge2,dataset2')
    result1,label_line1,second_line1=read_label_file.read(clf1)
    result2,label_line2,second_line2=read_label_file.read(clf2)
    second_line1.extend(second_line2)
    label=[str(int(i)+len(result1)) for i in label_line2]
    label_line1.extend(label)
    read_label_file.write(outfile,second_line1,label_line1)
    module_utils.write_Betsy_parameters_file(parameters,[obj1,obj2],pipeline)
    return new_objects

def make_unique_hash(parameters,objects):
    return module_utils.make_unique_hash(parameters,objects,'class_label_file','merge1,dataset1')
    
def get_outfile(parameters,objects):
    outfile = os.path.join(os.getcwd(),'class_label_file.cls')
    attributes = parameters.values()
    new_object=rule_engine.DataObject('class_label_file',attributes,outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return outfile,new_objects


