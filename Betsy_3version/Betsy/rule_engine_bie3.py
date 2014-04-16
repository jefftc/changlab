#rule_engine_bie3.py
import rulebase
import config
import bie3
import logging
import sys
import tempfile
import shutil
import module_utils
import os
import time
import getpass
from time import strftime,localtime

class DataObject:
    def __init__(self, data, identifier=""):
        self.data = data
        self.identifier = identifier
    def __repr__(self):
        x = str(self.data) + ' identifier:' + self.identifier
        return x
    
def compare_two_dict(dict_A,dict_B):
    dict_A_copy = dict_A.copy()
    dict_B_copy = dict_B.copy()
    if 'sf_processing_step' in dict_A_copy:
        del dict_A_copy['sf_processing_step']
    if 'sf_processing_step' in dict_B_copy:
        del dict_B_copy['sf_processing_step']
    if 'psf_processing_step' in dict_A_copy:
        del dict_A_copy['psf_processing_step']
    if 'psf_processing_step' in dict_B_copy:
        del dict_B_copy['psf_processing_step']
    if len(dict_A_copy)!=len(dict_B_copy):
        return False
    if set(dict_A_copy)-set(dict_B_copy):
        return False
    for key in dict_A_copy:
        setA = dict_A_copy[key]
        setB = dict_B_copy[key]
        if isinstance(setA,str):
            setA = [setA]
        if isinstance(setB,str):
            setB = [setB]
        if (not set(setA).issubset(set(setB))
            and not set(setB).issubset(set(setA))):
            return False
    return True
def get_outnode(network, module_id, module_file,parameters, pool,user_input):
    data_node = module_file.find_antecedents(network, module_id, pool)
    outfile = module_file.name_outfile(data_node,user_input)
    next_possible_ids = network.transitions[module_id]
    out_object = None
    parameters = module_file.get_out_attributes(parameters,data_node)
    fn=getattr(rulebase,network.nodes[next_possible_ids[0]].datatype.name)
    out_node = bie3.Data(fn,**parameters)
    out_object = DataObject(out_node,outfile)
    out_id = None
    for next_id in next_possible_ids:
        if compare_two_dict(parameters,network.nodes[next_id].attributes):
            out_id = next_id
##        next_node = create_object_with_full_attribute(network,
##                network.nodes[next_id], parameters,next_id)
##        current_attributes = parameters.copy()
##        new_attributes = next_node.attributes.copy()
##        keys = new_attributes.keys()
##        for key in keys:
##            DATA1_ATTR = network.nodes[next_id].datatype.get_attribute_object(key)
##            OPTIONAL = DATA1_ATTR.OPTIONAL
##            if OPTIONAL:
##                del current_attributes[key]
##                del new_attributes[key]
##        if cmp(new_attributes, current_attributes) == 0:
##            next_node.attributes = parameters
##            out_nodes.append((next_node, next_id))
        # e.g. 'v3' in ['v3',v4']
##        elif new_attributes.keys().sort() == current_attributes.keys().sort():
##            flag = True
##            for key in new_attributes:
##                if current_attributes[key] == new_attributes[key]:
##                    continue
##                elif (isinstance(new_attributes[key],list) and
##                      current_attributes[key] in new_attributes[key]):
##                    next_node.attributes[key] = current_attributes[key]
##                elif (isinstance(current_attributes[key],list) and
##                      new_attributes[key] in current_attributes[key]):
##                    next_node.attributes[key] = new_attributes[key]
##                elif current_attributes[key] == '___BETSY_ANYATOM___':
##                    next_node.attributes[key] = new_attributes[key]
##                else:
##                    flag = False
##            if flag:
##                out_nodes.append((next_node, next_id))
    print (out_object,out_id)
    return [(out_object,out_id)]

def create_out_attribtues_from_objects(network, module, module_id, pool):
    next_ids = network.transitions[module_id]
    #if len(next_ids)==1:
    return network.nodes[next_ids[0]].attributes
    #else:
        
    #return True
##    new_attributes = module.cons_data.attributes.copy()
##    for x in pool.keys():
##        node, node_id = pool[x],x
##        if module_id in network.transitions[node_id]:
##            attributes = node.attributes
##            defaults = module.cons_data.datatype.get_defaults()
##            #fill up the defaults parameters
##            for key in defaults.keys():
##                if key not in new_attributes:
##                    new_attributes[key]=defaults[key]
##            # pass the parameter from the previous Data node        
##            for key in attributes:
##                if key in defaults.keys():
##                    if new_attributes[key]==defaults[key]:
##                        new_attributes[key]=attributes[key]
##            # for module cons_data has attributes as a list but default is
##            #  used ,e.g. .preprocess_illumina
##            for key in defaults.keys():
##                if (isinstance(new_attributes[key],list) or
##                        new_attributes[key]=='___BETSY_ANYATOM___'):
##                        if defaults[key] in new_attributes[key]:
##                            new_attributes[key] = defaults[key]
##            break
##    # for values that can only be determine in the output node, like filter=25,
##    #['mean','median'] is shown
##    # in the outnode but not in the modules and previous nodes
##    out_ids = network.transitions[module_id]
##    for key in new_attributes:
##        if (new_attributes[key] == '___BETSY_ANYATOM___' or
##            isinstance(new_attributes[key],list)):
##            for out_id in out_ids:
##                out_node = network.nodes[out_id]
##                if out_node.attributes[key] != '___BETSY_ANYATOM___':
##                    new_attributes[key] = out_node.attributes[key]
##                    break
##    return new_attributes

def _can_module_take_data(module, datas):
    # Return True/False if a module can take this list of Data nodes
    # as an input.
    if len(module.in_datatypes) != len(datas):
        return False
    if len(datas) > 1:
        raise NotImplementedError
    data = datas[0]
    if data.datatype != module.in_datatypes[0]:
        return False

    # Make sure the data satisfies each of the module's constraints.
    for constraint in module.constraints:
        assert constraint.name in data.attributes
        data_value = data.attributes.get(constraint.name)
        data_type = bie3._get_attribute_type(data_value)
        assert data_type in [TYPE_ATOM, TYPE_ENUM]
        
        if constraint.behavior == MUST_BE:
            if data_type == TYPE_ATOM:
                if data_value != constraint.arg1:
                    return False
            elif data_type == TYPE_ENUM:
                return False
            else:
                raise AssertionError
        elif constraint.behavior == CAN_BE_ANY_OF:
            if data_type == TYPE_ATOM:
                if data_value not in constraint.arg1:
                    return False
            elif data_type == TYPE_ENUM:
                # data_value contains the possible values of this Data
                # object.  The values that are acceptable by module is
                # in constraint.arg1.  Make sure the module can handle
                # all of the possible values.
                if not bie3._is_subset(data_value, constraint.arg1):
                    return False
            else:
                raise AssertionError
        else:
            raise AssertionError

    return True


            
def _can_module_take_one_data_index(module, data, constraint_index):
    # Return True/False if a module can take this Data node
    # as part of the input with constraint_index.
        
    # Make sure the data satisfies each of the module's constraints.
    for constraint in module.constraints:
        #print constraint
        #ignore the sf_processing_step and psf_processing_step information
        if constraint.name in ['sf_processing_step','psf_processing_step']:
            continue
        if constraint.input_index != constraint_index:
            continue
        assert constraint.name in data.attributes
        data_value = data.attributes.get(constraint.name)
        data_type = bie3._get_attribute_type(data_value)
        assert data_type in [bie3.TYPE_ATOM, bie3.TYPE_ENUM]
        
        if constraint.behavior == bie3.MUST_BE:
            if data_type == bie3.TYPE_ATOM:
                if data_value != constraint.arg1:
                    return False
            elif data_type == bie3.TYPE_ENUM:
                return False
            else:
                raise AssertionError
        elif constraint.behavior == bie3.CAN_BE_ANY_OF:
            if data_type == bie3.TYPE_ATOM:
                if data_value not in constraint.arg1:
                    return False
            elif data_type == bie3.TYPE_ENUM:
                # data_value contains the possible values of this Data
                # object.  The values that are acceptable by module is
                # in constraint.arg1.  Make sure the module can handle
                # all of the possible values.
                if not bie3._is_subset(data_value, constraint.arg1):
                    return False
            else:
                raise AssertionError
        else:
            raise AssertionError

    return True
def _can_module_take_one_data(module, data):
    # Return True/False if a module can take this Data node
    # as part of the input.
    if data.datatype not in module.in_datatypes:
        return False
    constraint_index = 0
    while constraint_index < len(module.in_datatypes):
        flag = _can_module_take_one_data_index(module, data, constraint_index)
        if flag:
            return flag
        constraint_index = constraint_index+1
    return False

# choose the next module and return the id
def choose_next_module(network, node_id, node_data):
    if node_id not in network.transitions:
        return False
    next_node_ids = network.transitions[node_id]
    result = []
    for next_node_id in next_node_ids:
        next_node = network.nodes[next_node_id]
        assert isinstance(next_node,bie3.Module),'next node supposed to be a module'
        if _can_module_take_one_data(next_node,network.nodes[node_id]):
            result.append((next_node, next_node_id))
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
    require_data = []
    for node_id in require_id:
        
        if node_id not in pool:
            continue
        
        flag = _can_module_take_one_data(network.nodes[module_id],
                                         pool[node_id].data)
        if flag and node_id not in has_id:
            has_id.append(node_id)
    require_id.sort()
    has_id.sort()
    if require_id == has_id:
         return True
    common_id = list(set(require_id).intersection(set(has_id)))
    required_datatype = [data.name for data
                         in network.nodes[module_id].in_datatypes]
    common_datatype = [network.nodes[comid].datatype.name
                       for comid in common_id]
    required_datatype.sort()
    common_datatype.sort()
    if common_datatype == required_datatype:
        return True
    return False

def make_module_wd_name(network, module_file,module_id, 
                        user_input,pipeline_sequence,
                        pool,current_attributes):
    module_name = network.nodes[module_id].name
    data_object = module_file.find_antecedents(network, module_id, pool)
    hash_string = module_file.make_unique_hash(
            data_object, pipeline_sequence,current_attributes,user_input)
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
        
def run_module(network, module_id, pool, user_inputs, pipeline_sequence,
               user=getpass.getuser(), job_name='', clean_up=True):
    output_path = config.OUTPUTPATH
    assert os.path.exists(output_path), (
            'the output_path %s does not exist' % output_path)
    # get module
    module_node = network.nodes[module_id]
    sub_user_input = {}
    if module_node.user_inputs:
        for user_in in module_node.user_inputs:
            if user_in.name in user_inputs:
                sub_user_input[user_in.name] = user_inputs[user_in.name]
            if user_in.default is None:
                assert user_in.name in user_inputs,'%s should be specified' %user_in.name
    module_name = module_node.name
    print '[' + time.strftime('%l:%M%p') + '].' + module_name
    module = __import__('modules.' + module_name, globals(),
                        locals(), [module_name], -1)
    current_attributes = create_out_attribtues_from_objects(
        network, module_node, module_id, pool)
    # make directory name
    working_dir = os.path.join(output_path, make_module_wd_name(
        network, module, module_id, sub_user_input,pipeline_sequence,
        pool,current_attributes))
    # make name of outfile
    data_node = module.find_antecedents(network, module_id, pool)
    outfile = module.name_outfile(data_node,sub_user_input)
    temp_dir = ''
    temp_outfile = ''
    if not os.path.exists(os.path.join(working_dir,'stdout.txt')):
        #if no result has been generated, create temp folder and run analysis
        temp_dir = tempfile.mkdtemp()
        os.chdir(temp_dir)
        try:
            starttime = strftime(module_utils.FMT, localtime())
            out_node = module.run(data_node, current_attributes, sub_user_input, network)
            module_utils.write_Betsy_parameters_file(out_node.data.attributes,
                                                     data_node,sub_user_input,pipeline_sequence,
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
        network, module_id, module, current_attributes, pool,sub_user_input)
    assert out_nodes,'module %s fails' %module_node.name
    return out_nodes

def run_pipeline(network, in_objects, user_inputs, user=getpass.getuser(), job_name=''):
    output_path = config.OUTPUTPATH
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    LOG_FILENAME = os.path.join(output_path, 'traceback.txt')
    logging.basicConfig(filename=LOG_FILENAME,level=logging.DEBUG)
    try:
        pipeline_sequence = []
        stack_list = []
        pool = {}
        for data_object in in_objects:
            data_node = data_object.data
            start_node_ids = bie3._find_start_nodes(network,data_node)
            print 'start_node_ids',start_node_ids
            for start_node_id in start_node_ids:
                stack_list.append((data_object,start_node_id))
                pool[start_node_id]=data_object
        num_failures = 0
        print 'pool',pool
        while stack_list:
            assert num_failures < len(stack_list)
            stack_before_pop = stack_list[:]
            data_object, node_id = stack_list.pop()
            data_node = data_object
            if isinstance(data_node,DataObject):
                data_node = data_object.data
            if isinstance(data_node, bie3.Data):
                pool[node_id]=data_object
                result = choose_next_module(network, node_id, data_node)
                if result:
                   result.sort(key=lambda x:x[1])
                   for x in result:
                       module, module_id = x
                       stack_list.append((module, module_id))
            elif isinstance(data_node, bie3.Module):
                module_id = node_id
                test_required = test_require_data(network, module_id, pool)
                print 'module_id',module_id
                print 'test_required',test_required
                if test_required:       
                    pipeline_sequence.append(module.name)
                    out_nodes = run_module(network, module_id, pool, user_inputs, pipeline_sequence,
                                           user,job_name) ###
                    for x in out_nodes:
                        next_node, next_id = x
                        if next_id == 0:
                            if module_utils.exists_nz(next_node.identifier):
                                print ('['+ time.strftime('%l:%M%p') +
                                       '] Completed successfully and generated a file:')
                                print  next_node.identifier + '\r'
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
    except Exception, x:
            raise
    return True


