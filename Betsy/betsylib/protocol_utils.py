#protocol_utils.py
import os
import shutil
import Betsy_config
import hash_method

def get_result_folder(outfiles,parameters,pipeline):
    OUTPUTPATH = Betsy_config.OUTPUTPATH
    filename = os.path.split(outfiles[0][0])[-1]
    if '_BETSYHASH1_' in filename: 
        inputid = '_'.join(filename.split('_')[:-2])
    else:
        inputid = filename
    for i in range(len(outfiles[0])):
        folder_string = hash_method.hash_parameters(
            inputid,pipeline[0][i],**parameters[0][i])
        folder_name = 'result_folder_BETSYHASH1_'+folder_string
        result_folder = os.path.join(OUTPUTPATH,folder_name)
        if not os.path.exists(result_folder):
            os.mkdir(result_folder)
        for j in range(len(outfiles)):
            if len(outfiles[j]) == 1:
                final_output = os.path.split(outfiles[j][0])[-1]
                shutil.copyfile(outfiles[j][0],
                            os.path.join(result_folder,final_output))
            elif len(outfiles[j]) > 1:
                final_output = os.path.split(outfiles[j][i])[-1]
                shutil.copyfile(outfiles[j][i],
                            os.path.join(result_folder,final_output))
        f = file(os.path.join(result_folder,'Betsy_parameters.txt'),'w')
        f.write('Inputid:')
        f.write(inputid + '\n')
        f.write('Module output:\n')
        for key in parameters[0][i].keys():
             f.write(key+':'+parameters[0][i][key]+'\n')
        f.write('Pipeline module sequence:\n')
        f.write('\t'.join(pipeline[0][i]))
        f.close()

def format_prolog_query(
    predicate,dataset_id,content,parameters,modules):
    str_parameters = ','.join(parameters)
    output = str('['+dataset_id+'],[' + content + '],[' +
                 str_parameters + '],' + modules)
    query = predicate + '(' + output+')'
    return query

    
