#!/usr/bin/env python
#admin_utils.py
import os
import json
from stat import * 
import time
from Betsy import config,module_utils
import argparse
import shutil
from datetime import datetime
from time import strftime,localtime


output_path = config.OUTPUTPATH
depends_on_dict = dict()
generates_dict = dict()

def convert_unicode_str(text):
    """convert a unicode string to string(u'DatasetId' to 'DatasetId')"""
    if isinstance(text, list):
        final = []
        for subtext in text:
            final.append(convert_unicode_str(subtext))
    else:
        final = str(text)
    return final

def get_hash_key(filename):
    if '_BETSYHASH1_' not in filename:
        return filename
    filename = convert_unicode_str(filename)
    folder_name = os.path.split(os.path.split(filename)[-2])[-1]
    hash_string = folder_name.split('_')[-1]
    return hash_string

def get_fullname_from_key(key):
    result_folders = os.listdir(output_path)
    for result_folder in result_folders:
        if key in result_folder:
            return result_folder
    return None


def del_module_dependence(key,delete_list):
    global depends_on_dict
    global generates_dict
    if not delete_list:
        delete_list = [key]
    depends_list = depends_on_dict[key]
    for i in depends_list:
        generates = generates_dict[i]
        if generates==[key] and generates[0] not in delete_list:
            delete_list.append(generates[0])
            delete_list = del_module_dependence(generates[0],delete_list)
    return delete_list

            
def del_module_generates(key,delete_list):
    global depends_on_dict
    global generates_dict
    if not delete_list:
        delete_list = [key]
    if key in generates_dict:
        generates_list = generates_dict[key]
        for i in generates_list:
            if i not in delete_list:
                delete_list.append(i)
                delete_list = del_module_generates(i,delete_list)
    return delete_list



def get_size(start_path = '.'):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(start_path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            total_size += os.path.getsize(fp)
    return total_size

    
def get_status(key,depends_on_dict,generates_dict,starttime,modified_time,user,jobname):
    fullname = get_fullname_from_key(key)
    module_name = fullname.split('_BETSYHASH1')[0]
    dirname = os.path.join(output_path,fullname)
    st = os.stat(os.path.join(output_path,fullname))
    size = get_size(dirname)
    #modified_time = time.asctime(time.localtime(st[ST_MTIME]))
    FMT = module_utils.FMT
    a = datetime.strptime(modified_time,FMT)-datetime.strptime(starttime,FMT)
    run_time = str(a.days*24) +' h ' + str(a.seconds/60) +' m '+ str(a.seconds%60) +' s'
    dependence_hash = depends_on_dict[key]
    if key in generates_dict:
        generate_hash = generates_dict[key]
        result=[module_name,key,starttime,modified_time,run_time,size,dependence_hash,generate_hash,user,jobname]
    else:
        result=[module_name,key,starttime,modified_time,run_time,size,dependence_hash,[''],user,jobname]
    return result


def main():
    global depends_on_dict
    global generates_dict
    parser = argparse.ArgumentParser(description='run the protocol engine')
    parser.add_argument('--list_all',dest='list_all',action='store_const', default=False,
                        help='List all jobs', const=True)
    parser.add_argument('--del_depend', dest='del_depend', type=str,
                         default=None, help='del a job and dependences')
    parser.add_argument('--del_generate', dest='del_generate',
                        default=None,
                        type=str, help='del a job and everything it generates')
    parser.add_argument('--dry_run', dest='dry_run',default=False,
                        action='store_const',const=True,
                         help='show the folder that will be deleted')
    parser.add_argument('-u', dest='user',default=None,
                        type=str,
                         help='only manipulate the jobs by this user')
    parser.add_argument('--del_old', dest='del_old',default=30,
                        type=int,
                         help='delete the job n days old,default is 30 days')

    args = parser.parse_args()
    result_folders = os.listdir(output_path)
    depends_on_dict = dict()
    generates_dict = dict()
    final_module = []
    all_modules = []
    starttimes = []
    users = []
    jobnames=[]
    modified_times= []
    #generate depends_on_dict and generates_dict,all_modules,final_modules
    for result_folder in result_folders:
        if result_folder in ['._.DS_Store','.DS_Store','traceback.txt','userfile']:
            continue
        folder_path = os.path.join(output_path,result_folder)
        hash_string = result_folder.split('_')[-1]
        all_modules.append(hash_string)
        filenames = os.listdir(folder_path)
        assert 'Betsy_parameters.txt' in filenames,'no parameters file found in %s'%folder_path
        f = file(os.path.join(folder_path,'Betsy_parameters.txt'),'r')
        newline = f.read()
        f.close()
        jline = json.loads(newline)
        input_file = []
        starttimes.append(convert_unicode_str(jline[9]))
        modified_times.append(convert_unicode_str(jline[11]))
        users.append(convert_unicode_str(jline[13]))
        jobnames.append(convert_unicode_str(jline[15]))
        for j in jline[1]:
            if isinstance(j,list):
                input_file.extend([j[i*2+1] for i in range(len(j)/2)])
            else:
                input_file.extend([jline[1][i*2+1] for i in range(len(jline[1])/2)])
                break
        output_file = os.path.join(folder_path,jline[-1])
        if os.path.exists(output_file):
            key = get_hash_key(output_file)
            value = [get_hash_key(i) for i in input_file if os.path.exists(i)]
            depends_on_dict[key] = value
            for i in value:
                if i in generates_dict.keys():
                    generates_dict[i].append(key)
                else:
                    generates_dict[i] = [key]
    out_keys = depends_on_dict.keys()
    in_keys = generates_dict.keys()
    for i in out_keys:
        if i not in in_keys:
            final_module.append(i)
    result = []
    for i in range(len(all_modules)):
        module = all_modules[i]
        starttime = starttimes[i]
        user=users[i]
        jobname=jobnames[i]
        modified_time=modified_times[i]
        if args.user:
            if user == args.user:
                result.append(get_status(module,depends_on_dict,
                                         generates_dict,starttime,modified_time,user,jobname))
        else:
            result.append(get_status(module,depends_on_dict,
                                     generates_dict,starttime,modified_time,user,jobname))
    if args.list_all:
        column_name = ['Module name','Hash key','Start time','Finish time','Completion time',
                       'Space disk','Dependence key','Generates key','User','Job name']
        print '\t'.join(column_name)
        for i in result:
            i=[str(j) for j in i]
            print '\t'.join(i)
    if args.del_depend:
        del_key = args.del_depend
        assert del_key in all_modules,'cannot find the folder that contains the input string'
        del_list = del_module_dependence(del_key,[])
        folder_del = [get_fullname_from_key(i) for i in del_list]
        for i in folder_del:
            if args.dry_run:
                print os.path.join(output_path,i)
            else:
                shutil.rmtree(os.path.join(output_path,i))
    if args.del_generate:
        del_key = args.del_generate
        assert del_key in all_modules,'cannot find the folder that contains the input string'
        del_list = del_module_generates(del_key,[])
        folder_del = [get_fullname_from_key(i) for i in del_list]
        for i in folder_del:
            if args.dry_run:
                print os.path.join(output_path,i)
            else:
                shutil.rmtree(os.path.join(output_path,i))
    if args.del_old:
        FMT = module_utils.FMT
        current_time = strftime(FMT,localtime())
        for i in result:
            a = datetime.strptime(current_time,FMT)-datetime.strptime(i[2],FMT)
            if a.days > args.del_old:
                folder_name = get_fullname_from_key(i[1])
                shutil.rmtree(os.path.join(output_path,folder_name))

if __name__ == '__main__':
    main()
