from Betsy import rulebase
from Betsy import bie
import json
import itertools

def convert_to_builtin_type(obj):
    # Convert objects to a dictionary of their representation
    d = { '__class__':obj.__class__.__name__, 
          '__module__':obj.__module__,
          }
    d.update(obj.__dict__)
    return d


def dict_to_object(d):
    if '__class__' in d:
        inst=None
        class_name = d.pop('__class__')
        module_name = d.pop('__module__')
        if '.' in module_name:
            module = __import__(module_name, globals(),
                        locals(), [module_name.split('.')[-1]], -1)
        else:
            module = __import__(module_name)
        class_ = getattr(module, class_name)
        args = dict( (key.encode('ascii'), value) for key, value in d.items())
        for key in args:
            if isinstance(args[key],dict):
                args[key]=dict_to_object(args[key])
            if isinstance(args[key],list) and isinstance(args[key][0],unicode):
                args[key]=[i.encode('ascii') for i in args[key]]
            if isinstance(args[key],unicode):
                args[key]=args[key].encode('ascii')
        if class_name=='Attribute':
            if  'name' in args.keys():
                args[args['name'].encode('ascii')]=args['values']
            del args['name']
            del args['values']
            inst = class_(**args)
        elif class_name == 'DataType':
            inst = class_(args['name'],*args['attribute_objects'])
        elif class_name == 'Data':
            inst = class_(args['datatype'],**args['attributes'])
        elif class_name == 'Network':
            new_transition = dict()
            for key,value in args['transitions'].items():
                new_transition[int(key)]=value
            inst = class_(args['nodes'],new_transition)
        elif class_name in ['Module','QueryModule']:
            inst = class_(args['ante_datas'],args['cons_data'],args['name'])
        else:
            assert 'else module %s' %class_name
    else:
        args = dict( (key.encode('ascii'), value) for key, value in d.items())
        for key in args:
            if isinstance(args[key],dict):
                args[key]=dict_to_object(args[key])
            if isinstance(args[key],list) and isinstance(args[key][0],unicode):
                args[key]=[i.encode('ascii') for i in args[key]]
            if isinstance(args[key],unicode):
                args[key]=args[key].encode('ascii')
        inst = args
    return inst

def get_network_list():
	attributes =  rulebase.SignalFile2.attributes
	defaults= rulebase.SignalFile2.get_defaults()
	all_list=[]
	key_name=[]
	for key in attributes:
	    if key=='filename':
		continue
            if attributes[key]=='___BETSY_ANYATOM___':
        	attributes[key]=['ANY',defaults[key]]
   	    key_name.append(key)
            all_list.append(attributes[key])
       
	comb=itertools.product(*all_list)
	goal_datatype=rulebase.SignalFile2
	for j in range(10):
            print j
            final_list=[]
            for i in itertools.islice(comb,100*j,100*(j+1)):
                goal_attributes=dict()
                for index,key in enumerate(key_name):
                    goal_attributes[key]=i[index]
                parameter = goal_attributes.copy()
                network = bie.backchain(
                    rulebase.all_modules, goal_datatype, goal_attributes)
    ##    	    network = bie.select_start_node(network, in_data)
                network = bie.optimize_network(network)
                if len(network.nodes)>1:
                    final_list.append(network)
            filename = 'pipeline'+str(j)+'.txt'
            f = file(filename,'w')
            json.dump(final_list,f,default=convert_to_builtin_type,indent=2)
            f.flush()
            print 'done save'
            f = file(filename,'r')
            text =f.read()
            f.close()
            myobj_instance = json.loads(text, object_hook=dict_to_object)
            print 'len',len(myobj_instance)


get_network_list()

