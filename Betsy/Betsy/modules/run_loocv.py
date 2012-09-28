#run_loocv.py
import arrayio
from Betsy import read_label_file
from Betsy import module_utils
import os
import svmutil

def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    
    outfile = get_outfile(parameters,objects,pipeline)
    M = arrayio.read(single_object.identifier)
    training_label_file = module_utils.find_object(parameters,
                                    objects,'class_label_file','contents')
    assert os.path.exists(training_label_file.identifier),(
        'the training label file %s does not exist'
        %training_label_file.identifier)
    a,training_label,second_line = read_label_file.read(training_label_file.identifier)
    full_index = range(M.ncol())
    train_model = None
    if 'train_model' in parameters.keys():
        train_model = __import__(parameters['train_model'])
    predict_model = __import__(parameters['predict'])
    f = file(outfile,'w')
    f.write('\t'.join(['sample_name','Predicted_class','Confidence','Actual_class','Correct?']))
    f.write('\n')
    for i in range(M.ncol()):
        test_index = i
        train_index = full_index[:]
        train_index.remove(test_index)
        y_training = [training_label[x] for x in train_index]
        y_test = [training_label[test_index]]
        M_train = M.matrix(None,train_index)
        M_test = M.matrix(None,[test_index])
        train_file = 'train'+'_'+str(i)
        test_file = 'test'+'_'+str(i)
        if parameters['format'] == 'tdf':
            f_out = file(train_file,'w')
            arrayio.tab_delimited_format.write(M_train,f_out)
            f_out.close()
            f_out = file(test_file,'w')
            arrayio.tab_delimited_format.write(M_test,f_out)
            f_out.close()
        elif parameters['format'] == 'gct':
            f_out = file(train_file,'w')
            arrayio.gct_format.write(M_train,f_out)
            f_out.close()
            f_out = file(test_file,'w')
            arrayio.gct_format.write(M_test,f_out)
            f_out.close()
        train_label = 'train_label'+'_'+ str(i)
        test_label = 'test_label' + '_' + str(i)
        read_label_file.write(train_label,second_line,y_training)
        read_label_file.write(test_label,[second_line[int(y_test[0])]],'0')
        train_parameters = parameters.copy()
        train_parameters['contents'] = '[train]'
        train_attributes = train_parameters.values()
        test_parameters = parameters.copy()
        test_parameters['contents'] = '[test]'
        test_attributes = test_parameters.values()
        train_object = module_utils.DataObject('signal_file',train_attributes,train_file)
        test_object = module_utils.DataObject('signal_file',test_attributes,test_file)
        train_clfobject = module_utils.DataObject('class_label_file',train_attributes,train_label)
        test_attributes.append('given')
        test_clfobject = module_utils.DataObject('class_label_file',test_attributes,test_label)
        temp_objects = [train_object,test_object,train_clfobject,test_clfobject]
        loo_parameters = parameters.copy()
        if train_model:
            loo_parameters['contents'] = '[train]'
            temp_objects = train_model.run(loo_parameters,temp_objects,[])
        temp_file = temp_objects[-1].identifier
        loo_parameters['traincontents'] ='[train]'
        loo_parameters['testcontents'] ='[test]'
        out_objects = predict_model.run(loo_parameters,temp_objects,[])
        loo_outfile = out_objects[-1].identifier
        f1 = open(loo_outfile,'r')
        lines = f1.readlines()
        f1.close()
        f.write(lines[1])
        os.remove(train_file)
        os.remove(test_file)
        os.remove(train_label)
        os.remove(test_label)
        os.remove(loo_outfile)
        if train_model:
            os.remove(temp_file)
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for loocv fails'%outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline,outfile)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)
    
def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'predication_loocv_'+original_file+'.txt'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for loocv does not exist'
        %single_object.identifier)
    return single_object


def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,'classification_file',parameters,objects,single_object)
    return new_objects
