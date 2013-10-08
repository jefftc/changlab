#rank_genes_by_class_neighbors.py
#from Betsy
import module_utils
import shutil
import os
from genomicode import config
import subprocess
from time import strftime,localtime
import bie
import rulebase

def run(in_nodes, parameters):
    data_node,cls_node = in_nodes
    outfile = name_outfile(in_nodes)
    import arrayio
    tmp = os.path.join(os.getcwd(),'tmp.txt')
    f = file(tmp,'w')
    M = arrayio.read(data_node.attributes['filename'])
    M_c = arrayio.convert(M,to_format=arrayio.gct_format)
    arrayio.gct_format.write(M_c,f)
    f.close()
    module_name = 'ClassNeighbors'
    gp_parameters = dict()
    gp_parameters['data.filename'] = tmp
    gp_parameters['class.filename'] = cls_node.attributes['filename']
    gp_parameters['num.neighbors'] = str(parameters['cn_num_neighbors'])
    
    if  parameters['cn_num_perm'].isdigit():
        gp_parameters['num.permutations'] = str(parameters['cn_num_perm'])
    if  module_utils.is_number(parameters['cn_user_pval']):
        gp_parameters['user.pval'] = str(parameters['cn_user_pval'])
        
    mean_median = {'cn_mean':'','cn_median':'-d'}
    if parameters['cn_mean_or_median'] in ['cn_mean','cn_median']:
        gp_parameters['mean.or.median'] = mean_median[parameters['cn_mean_or_median']]
        
    p={'cn_ttest':'','cn_snr':'-S'}
    if  parameters['cn_ttest_or_snr'] in p.values():
        gp_parameters['ttest.or.snr'] = p[parameters['cn_ttest_or_snr']]
        
    if  parameters['cn_ttest_or_snr'] in ['cn_yes','cn_no']:
        gp_parameters['filter.data'] = str(parameters['cn_filter_data'])[3:]

    if  module_utils.is_number(parameters['cn_abs_diff']):
        gp_parameters['min.abs.diff'] = str(parameters['cn_abs_diff'])
        
    if  module_utils.is_number(parameters['cn_min_threshold']):
        gp_parameters['min.threshold'] = str(parameters['cn_min_threshold'])
        
    if  module_utils.is_number(parameters['cn_max_threshold']):
        gp_parameters['max.threshold'] = str(parameters['cn_max_threshold'])
        
    if  module_utils.is_number(parameters['cn_min_folddiff']):
        gp_parameters['min.fold.diff'] = str(parameters['cn_min_folddiff'])
        
    gp_path = config.genepattern
    gp_module = module_utils.which(gp_path)
    assert gp_module,'cannot find the %s' %gp_path
    download_directory = os.path.join(os.getcwd(),'class_neighbors_result')
    command = [gp_module, module_name,'-o', download_directory]
    for key in gp_parameters.keys():
        a = ['--parameters',key+':'+ gp_parameters[key]]
        command.extend(a)
    process = subprocess.Popen(command,shell=False,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    process.wait()
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
        'the output file %s for rank_genes_by_class_neighbors fails' %outfile)
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.GeneListFile,**new_parameters)
    return out_node

def find_antecedents(network, module_id,data_nodes):
    data_node = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='SignalFile2')
    cls_node = module_utils.get_identifier(network, module_id, data_nodes,
                                           datatype='ClassLabelFile')
    return data_node, cls_node

def name_outfile(in_nodes):
    data_node,cls_node = in_nodes
    original_file = module_utils.get_inputid(
        data_node.attributes['filename'])
    filename = 'gene_list' + original_file + '.txt'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters):
    data_node,cls_node = in_nodes
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)




