#rank_gene_by_class_neighbors.py
from Betsy import module_utils
import shutil
import os
from genomicode import jmath
from Betsy import config
import subprocess

def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    label_file = module_utils.find_object(
        parameters,objects,'class_label_file','contents')
    assert os.path.exists(label_file.identifier),(
        'cannot find label_file %s for class_neighbors'%label_file.identifier)
    import arrayio
    tmp = 'tmp.txt'
    f = file(tmp,'w')
    M = arrayio.read(single_object.identifier)
    M_c = arrayio.convert(M,to_format=arrayio.gct_format)
    arrayio.gct_format.write(M_c,f)
    f.close()
    module_name = 'ClassNeighbors'
    gp_parameters = dict()
    gp_parameters['data.filename'] = tmp
    gp_parameters['class.filename'] = label_file.identifier
    
    if 'cn_num_neighbors' in parameters.keys():
        assert parameters['cn_num_neighbors'].isdigit(),'cn_num_neighbors should be digit'
        gp_parameters['num.neighbors'] = str(parameters['cn_num_neighbors'])
        
    if 'cn_num_perm' in parameters.keys():
        assert  gp_parameters['cn_num_perm'].isdigit(),'cn_num_perm should be digit'
        gp_parameters['num.permutations'] = str(parameters['cn_num_perm'])
        
    if 'cn_user_pval' in parameters.keys():
        assert  module_utils.is_number(gp_parameters['cn_user_pval']),'cn_user_pval should be a number'
        gp_parameters['user.pval'] = str(parameters['cn_user_pval'])
        
    if 'cn_mean_or_median' in parameters.keys():
        mean_median = {'cn_mean':'','cn_median':'-d'}
        assert  gp_parameters['cn_mean_or_median'] in ['cn_mean','cn_median'],'ill_coll_mode is not correct'
        gp_parameters['mean.or.median'] = mean_median[parameters['cn_mean_or_median']]
        
    if 'cn_ttest_or_snr' in parameters.keys():
        p={'cn_ttest':'','cn_snr':'-S'}
        assert parameters['cn_ttest_or_snr'] in p.values(),'cn_ttest_snr is invalid'
        gp_parameters['ttest.or.snr'] = p[parameters['cn_ttest_or_snr']]
        
    if 'cn_filter_data' in parameters.keys():
        assert parameters['cn_ttest_or_snr'] in ['cn_yes','cn_no'],'cn_filter_data is invalid'
        gp_parameters['filter.data'] = str(parameters['cn_filter_data'])[3:]

    if 'cn_abs_diff' in parameters.keys():
        assert module_utils.isnumber(parameters['cn_abs_diff']),'cn_abs_diff should be number'
        gp_parameters['min.abs.diff'] = str(parameters['cn_abs_diff'])
        
    if 'cn_min_threshold' in parameters.keys():
        assert module_utils.isnumber(parameters['cn_min_threshold']),'cn_min_threshold should be number'
        gp_parameters['min.threshold'] = str(parameters['cn_min_threshold'])
        
    if 'cn_max_threshold' in parameters.keys():
        assert module_utils.isnumber(parameters['cn_max_threshold']),'cn_max_threshold should be number'
        gp_parameters['max.threshold'] = str(parameters['cn_max_threshold'])
        
    if 'cn_min_folddiff' in parameters.keys():
        assert module_utils.isnumber(parameters['cn_min_folddiff']),'cn_min_folddiff should be number'
        gp_parameters['min.fold.diff'] = str(parameters['cn_min_folddiff'])
        
    gp_path = config.GENEPATTERN
    gp_module = module_utils.which(gp_path)
    assert gp_module,'cannot find the %s' %gp_path
    command = [gp_module, module_name]
    for key in gp_parameters.keys():
        a = ['--parameters',key+':'+ gp_parameters[key]]
        command.extend(a)
    
    download_directory = None
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    process.wait()
    out_text =  process.stdout.read()
    out_lines = out_text.split('\n')
    for out_line in out_lines:
        if out_line != 'Loading required package: rJava' and len(out_line)>0:
            download_directory = out_line
            break
    error_message = process.communicate()[1]
    if error_message:
        raise ValueError(error_message)
    assert os.path.exists(download_directory),(
        'there is no output directory for class_neighbors')
    result_files = os.listdir(download_directory)
    assert 'stderr.txt' not in result_files,'gene_pattern get error'
    os.remove(tmp)
    gene_list=[]
    for result_file in result_files:
        if result_file.endswith('.odf'):
            f = file(os.path.join(download_directory,result_file),'r')
            text = f.read()
            text = text.split('\n')
            f.close()
            numline = 8
            startline = 14
            assert text[numline].startswith('NumNeighbors'),'the odf file format is not right'
            number_gene=int(text[numline].split('=')[1])
            assert text[startline].startswith('1'),'the start line is not right'
            
            for line in text[startline:startline+number_gene]:
                lines = line.split('\t')
                gene_list.append(lines[10])
    f = file(outfile,'w')
    f.write('\t'.join(gene_list))
    f.close()
    assert module_utils.exists_nz(outfile),(
        'the output file %s for class_neighbors fails' %outfile)
    new_objects = get_newobjects(parameters,objects,pipeline)
    module_utils.write_Betsy_parameters_file(
        parameters,single_object,pipeline,outfile)
    return new_objects

def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,'signal_file','contents,preprocess')
    assert os.path.exists(single_object.identifier),(
        'the input file %s for class_neighbors does not exist'
        %single_object.identifier)
    return single_object

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = 'gene_list_'+original_file+'.txt'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile
    
def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    parameters = module_utils.renew_parameters(parameters,['status'])
    new_object = module_utils.DataObject(
        'gene_list_file',[parameters['contents']],outfile)
    new_objects = objects[:]
    new_objects.append(new_object)
    return new_objects
