#rule_engine_bie.py
### rule engine with filename as an attribute in Data,

import config
import getpass
import logging
import sys
import tempfile
import shutil
import module_utils
import bie
import os
import time
from time import strftime,localtime
# Functions:
# make_module_wd_name
# make_module_wd
# copy_result_folder

# _is_antecedence_compatible_with_module
# choose_next_module
# test_require_data
# create_object_with_full_attribute
# create_module_attribtues_from_objects
# get_outnode

# run_module
# run_pipeline



def get_outnode(network, module_id, module_file,parameters, pool):
    data_node = module_file.find_antecedents(network, module_id, pool,parameters)
    outfile = module_file.name_outfile(data_node)
    next_possible_ids = network.transitions[module_id]
    parameters = module_file.get_out_attributes(parameters,data_node)
    parameters['filename'] = outfile
    out_nodes = []
    for next_id in next_possible_ids:
        next_node = create_object_with_full_attribute(network,
                network.nodes[next_id], parameters,next_id)
        current_attributes = parameters.copy()
        new_attributes = next_node.attributes.copy()
        keys = new_attributes.keys()
        for key in keys:
            DATA1_ATTR = network.nodes[next_id].datatype.get_attribute_object(key)
            OPTIONAL = DATA1_ATTR.OPTIONAL
            if OPTIONAL:
                del current_attributes[key]
                del new_attributes[key]
        if cmp(new_attributes, current_attributes) == 0:
            next_node.attributes = parameters
            out_nodes.append((next_node, next_id))
        # e.g. 'v3' in ['v3',v4']
        elif new_attributes.keys().sort() == current_attributes.keys().sort():
            flag = True
            for key in new_attributes:
                if current_attributes[key] == new_attributes[key]:
                    continue
                elif (isinstance(new_attributes[key],list) and
                      current_attributes[key] in new_attributes[key]):
                    next_node.attributes[key] = current_attributes[key]
                elif (isinstance(current_attributes[key],list) and
                      new_attributes[key] in current_attributes[key]):
                    next_node.attributes[key] = new_attributes[key]
                elif current_attributes[key] == '___BETSY_ANYATOM___':
                    next_node.attributes[key] = new_attributes[key]
                else:
                    flag = False
            if flag:
                out_nodes.append((next_node, next_id))
    return out_nodes

def make_module_wd_name(network, module_file, module_id, current_attributes,
                        pipeline_sequence, pool):
    module_name = network.nodes[module_id].name
    data_node = module_file.find_antecedents(
        network, module_id, pool,current_attributes)
    hash_string = module_file.make_unique_hash(
            data_node, pipeline_sequence,
            current_attributes)
    working_dir = module_name + '_BETSYHASH1_' + hash_string
    return working_dir


def make_module_wd(working_dir):
    os.mkdir(working_dir)

def copy_result_folder(working_dir, temp_dir, temp_outfile):
    """copy result files in temp folder to result folder,
      if fails, delete result folder"""
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)
    try:
        shutil.copyfile(os.path.join(temp_dir,'Betsy_parameters.txt'),
                        os.path.join(working_dir,'Betsy_parameters.txt'))
        if os.path.isdir(temp_outfile):
            shutil.copytree(os.path.join(temp_dir,temp_outfile),
                                     os.path.join(working_dir,temp_outfile))
        else:
            #change to copy other files
            all_files = os.listdir(temp_dir)
            extra_files = list(set(all_files)-set([temp_outfile,
                                                   'Betsy_parameters.txt','stdout.txt']))
            shutil.copyfile(os.path.join(temp_dir,temp_outfile),
                                     os.path.join(working_dir,temp_outfile))
            for extra_file in extra_files:
                if os.path.isfile(os.path.join(temp_dir,extra_file)):
                    shutil.copyfile(os.path.join(temp_dir,extra_file),
                                     os.path.join(working_dir,extra_file))
                elif os.path.isdir(os.path.join(temp_dir,extra_file)):
                    shutil.copytree(os.path.join(temp_dir,extra_file),
                                     os.path.join(working_dir,extra_file))
        shutil.copyfile(os.path.join(temp_dir,'stdout.txt'),
                                     os.path.join(working_dir,'stdout.txt'))
    finally:
        if not os.path.isfile(os.path.join(working_dir,'stdout.txt')):
            shutil.rmtree(working_dir)
            raise ValueError('copy fails')
        
def is_antecedence_compatible_with_module(data, module):
    # data is a Data node in the network. module is the data
    # connect in the network
    ante_datatype = [x.datatype for x in module.ante_datas]
    if data.datatype not in ante_datatype:
        return False
    index = 0
    ante_datatype = [x for x in ante_datatype if x == data.datatype]
    temp_ante_datas = [x for x in module.ante_datas if x.datatype == data.datatype]
    #consider merge multiple SignalFile and SignatureScore
    if len(ante_datatype)>1 and ante_datatype == [ante_datatype[0]]*len(ante_datatype):
        data_keys = data.attributes.keys()
        for i, x in enumerate(temp_ante_datas):
            common_keys = set(x.attributes.keys()).intersection(set(data_keys))
            flag = True
            for common_key in common_keys:
                if x.attributes[common_key] != data.attributes[common_key]:
                    flag = False
                    break
            if  flag:
                index = i
                break
    else:
        index = ante_datatype.index(data.datatype)
    data_attr = data.attributes
    module_attr = temp_ante_datas[index].attributes
    compatible = True
    # module antecendent Data
    # ATOM      Must match a specific value.
    # ENUM      UNDEFINED.
    # ANYATOM   Can match any value.
    # NOVALUE   Use default value.
    # 
    # CASE ANTE_TYPE  DATA_TYPE   RESULT
    #   1    NOVALUE    NOVALUE    no.
    #   2    NOVALUE    ANYATOM    no.  
    #   3    NOVALUE      ATOM     no.
    #   4    NOVALUE      ENUM     no
    #   5    ANYATOM    NOVALUE    ERROR.  DATA should have default values.
    #   6    ANYATOM    ANYATOM    OK.
    #   7    ANYATOM      ATOM     OK.
    #   8    ANYATOM      ENUM    No.
    #   9      ATOM     NOVALUE    ERROR.  DATA should have default values
    #  10      ATOM     ANYATOM    OK.
    #  11      ATOM       ATOM     Check if items are equal.
    #  12      ATOM       ENUM     Check if ATOM in ENUM.
    #  13      ENUM     NOVALUE    ERROR.  DATA should have default values.
    #  14      ENUM     ANYATOM    NotImplementedError.
    #  15      ENUM       ATOM     Check if ATOM in ENUM.
    #  16      ENUM       ENUM     check if some intersection.
    for key in module_attr:
        DATA_VALUE = data_attr.get(key)
        ANTE_VALUE = module_attr.get(key)
        DATA_TYPE = bie._get_attribute_type(data_attr, key)
        ANTE_TYPE = bie._get_attribute_type(module_attr, key)
        CASE = bie._assign_case_by_type(ANTE_TYPE, DATA_TYPE)
        if CASE in [5, 9, 13]:
            raise AssertionError
        elif CASE in [14]:
            raise NotImplementedError
        elif CASE in [1, 2, 3, 4, 8]:
            compatible = False
        elif CASE in [6, 7, 10]:
            pass
        elif CASE == 11:
            if DATA_VALUE != ANTE_VALUE:
                compatible = False
        elif CASE == 12:
            if ANTE_VALUE not in DATA_VALUE:
                compatible = False
        elif CASE == 15:
            if DATA_VALUE not in ANTE_VALUE:
                compatible = False
        elif CASE == 16:
            if not set(DATA_VALUE).intersection(set(ANTE_VALUE)):
                compatible = False
        else:
            raise AssertionError
    return compatible
    


# choose the next module and return the id
def choose_next_module(network, node_id, node_data):
    if node_id not in network.transitions:
        return False
    next_node_ids = network.transitions[node_id]
    result = []
    for next_node_id in next_node_ids:
        next_node = network.nodes[next_node_id]
        assert isinstance(next_node,bie.Module),'next node supposed to be a module'
        module_cons_node = create_object_with_full_attribute(network,
            next_node.cons_data, node_data.attributes)
        new_module = bie.Module(next_node.name, next_node.ante_datas[:], module_cons_node)
        if is_antecedence_compatible_with_module(network.nodes[node_id], next_node):
            result.append((new_module, next_node_id))
    result.sort(key=lambda x: x[1], reverse=False)
    return result



def test_require_data(network, module_id, pool):
    # test if the required data for the module is met.
    require_id = []
    data_type = []
    has_id = []
    for key in network.transitions:
        if module_id in network.transitions[key]:
            require_id.append(key)
    for node_id in require_id:
        if node_id not in pool:
            continue
        flag = bie._is_compatible_with_start(network.nodes[node_id],pool[node_id])
        if flag and node_id not in has_id:
            has_id.append(node_id)
    require_id.sort()
    has_id.sort()
    if require_id == has_id:
         return True
    common_id = list(set(require_id).intersection(set(has_id)))
    required_datatype = [data.datatype.name for data
                         in network.nodes[module_id].ante_datas]
    common_datatype = [network.nodes[comid].datatype.name for comid in common_id]
    required_datatype.sort()
    common_datatype.sort()
    if common_datatype == required_datatype:
        return True
    return False



def create_object_with_full_attribute(network,in_data,ante_attributes,node_id=None):
    new_data = bie.Data(in_data.datatype, **in_data.attributes)
    attributes = new_data.attributes
    if node_id:
        attributes_in_network = network.nodes[node_id].attributes
        # e.g. pass GSEID to GEOSeries
        for key in attributes_in_network:
            if key not in attributes:
                new_data.attributes[key] = attributes_in_network[key]
    for key in ante_attributes:
        if key in in_data.datatype.get_defaults():
            if key not in attributes:
                new_data.attributes[key] = ante_attributes[key]
            elif (isinstance(attributes[key],list) and
                  ante_attributes[key] in attributes[key]):
                new_data.attributes[key] = ante_attributes[key]
    return new_data



def create_module_attribtues_from_objects(network, module, module_id, pool):
    new_attributes = module.cons_data.attributes.copy()
    for x in pool.keys():
        node, node_id = pool[x],x
        if module_id in network.transitions[node_id]:
            attributes = node.attributes
            defaults = module.cons_data.datatype.get_defaults()
            #fill up the defaults parameters
            for key in defaults.keys():
                if key not in new_attributes:
                    new_attributes[key]=defaults[key]
            # pass the parameter from the previous Data node        
            for key in attributes:
                if key in defaults.keys():
                    if new_attributes[key]==defaults[key]:
                        new_attributes[key]=attributes[key]
            # for module cons_data has attributes as a list but default is
            #  used ,e.g. .preprocess_illumina
            for key in defaults.keys():
                if (isinstance(new_attributes[key],list) or
                        new_attributes[key]=='___BETSY_ANYATOM___'):
                        if defaults[key] in new_attributes[key]:
                            new_attributes[key] = defaults[key]
            break
    # for values that can only be determine in the output node, like filter=25,
    #['mean','median'] is shown
    # in the outnode but not in the modules and previous nodes
    out_ids = network.transitions[module_id]
    for key in new_attributes:
        if (new_attributes[key] == '___BETSY_ANYATOM___' or
            isinstance(new_attributes[key],list)):
            for out_id in out_ids:
                out_node = network.nodes[out_id]
                if out_node.attributes[key] != '___BETSY_ANYATOM___':
                    new_attributes[key] = out_node.attributes[key]
                    break
    return new_attributes


def run_module(network, module_id, module_node, pool, pipeline_sequence,
               user=getpass.getuser(), job_name='', clean_up=True):
    output_path = config.OUTPUTPATH
    assert os.path.exists(output_path), (
            'the output_path %s does not exist' % output_path)
    # get module
    module_name = module_node.name
    print '[' + time.strftime('%l:%M%p') + '].' + module_name
    module = __import__('modules.' + module_name, globals(),
                        locals(), [module_name], -1)
    current_attributes = create_module_attribtues_from_objects(
        network, module_node, module_id, pool)
    # make directory name
    working_dir = os.path.join(output_path, make_module_wd_name(
        network, module, module_id, current_attributes,
        pipeline_sequence, pool))
    # make name of outfile
    data_node = module.find_antecedents(network, module_id, pool,current_attributes)
    outfile = module.name_outfile(data_node)
    temp_dir = ''
    temp_outfile = ''
    if not os.path.exists(os.path.join(working_dir,'stdout.txt')):
        #if no result has been generated, create temp folder and run analysis
        temp_dir = tempfile.mkdtemp()
        os.chdir(temp_dir)
        try:
            starttime = strftime(module_utils.FMT, localtime())
            out_node = module.run(data_node, current_attributes, network)
            module_utils.write_Betsy_parameters_file(out_node.attributes,
                                                     data_node,pipeline_sequence,
                                                     starttime,user,job_name)
            if out_node:
                f = file(os.path.join(temp_dir, 'stdout.txt'),'w')
                f.write('This module runs successully.')
                f.close()
            else:
                return False
        except:
            logging.exception('Got exception on %s .run()'%module_name)
            raise
        temp_outfile = os.path.split(outfile)[-1]
    try:
        # found if the same analysis has been run and wait for it finished
        try:
            os.mkdir(working_dir)
        except OSError:
            try :
                module_timeout = int(config.MODULE_TIMEOUT)
            except AttributeError:
                module_timeout = 60
            i = 0
            while i <= module_timeout:
                if os.path.isfile(os.path.join(working_dir,'stdout.txt')):
                    break
                i = i + 1
                time.sleep(1)
        # if previous analysis not working, copy the current results to working_dir
        if (module_utils.exists_nz(temp_outfile) and not
            os.path.exists(os.path.join(working_dir,'stdout.txt'))):
            copy_result_folder(working_dir, temp_dir, temp_outfile)
    finally:
        if clean_up and os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)
    os.chdir(working_dir)
    out_nodes = get_outnode(
        network, module_id, module, current_attributes, pool)
    assert out_nodes,'module %s fails' %module_node.name
    return out_nodes

def run_pipeline(network, in_data, user=getpass.getuser(), job_name=''):
    output_path = config.OUTPUTPATH
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    LOG_FILENAME = os.path.join(output_path, 'traceback.txt')
    logging.basicConfig(filename=LOG_FILENAME,level=logging.DEBUG)
    try:
        pipeline_sequence = []
        stack_list = []
        pool = {}
        for data_node in in_data:
            start_node_ids = bie._find_start_nodes(network, data_node)
            for start_node_id in start_node_ids:
                new_data = create_object_with_full_attribute(
                    network,data_node, data_node.datatype.get_defaults(),start_node_id)
                stack_list.append((new_data,start_node_id))
        num_failures = 0
        while stack_list:
            assert num_failures < len(stack_list)
            stack_before_pop = stack_list[:]
            data_node, node_id = stack_list.pop()
            if isinstance(data_node, bie.Data):
                pool[node_id]=data_node
                result = choose_next_module(network, node_id, data_node)
                if result:
                   result.sort(key=lambda x:x[1])
                   for x in result:
                       module, module_id = x
                       stack_list.append((module, module_id))
            elif isinstance(data_node, bie.Module):
                module = data_node
                module_id = node_id
                test_required = test_require_data(network, module_id, pool)
                if test_required:       
                    pipeline_sequence.append(module.name)
                    out_nodes = run_module(network, module_id, module, pool, pipeline_sequence,
                                           user,job_name)
                    for x in out_nodes:
                        next_node, next_id = x
                        if next_id == 0:
                            if module_utils.exists_nz(next_node.attributes['filename']):
                                print ('['+ time.strftime('%l:%M%p') +
                                       '] Completed successfully and generated a file:')
                                print  next_node.attributes['filename'] + '\r'
                                print '\r'
                                sys.stdout.flush()
                            else:
                                    print 'This pipeline has completed unsuccessfully'
                                    raise ValueError(
                                        'there is no output for this pipeline')
                            return
                        stack_list.append((next_node,next_id))
                else:
                    if stack_list:
                        stack_list.insert(0,(data_node, node_id))
                        num_failures += 1
        if not stack_list and node_id > 0 : # if stack_list empty but have not find node 0
            assert ValueError('cannot find the final node')
    except Exception, x:
            raise
    return True




