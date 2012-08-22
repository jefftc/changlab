#arrayannot.py

import os
import filelib,arrayplatformlib
import config
import re
from arrayio import const
import arrayio

def hash_chipname(filename):
    x = os.path.split(filename)[1]
    x = x.replace(".gz", "")
    x = x.replace(".csv", "")
    x = x.replace(".annot", "")
    x = x.replace("_annot", "")
    x = x.replace("-", "_")
    version = 0
    m = re.search(r".na(\d+)",x)
    if m:
        version = m.group(1)
        x = x.replace(m.group(0),'')
    return x,version

def chipname2filename(chipname):
    filename = chipname2filename_affy(chipname)
    if not filename:
        filename = chipname2filename_illu(chipname)
    return filename

def chipname2filename_illu(chipname):
    filename = None
    path = config.annot_data_iilu
    assert os.path.exists(path),'%s does not exist'%path
    chipname=chipname.replace('_','-')
    for file in os.listdir(path):
        if chipname in file:
            filename = os.path.join(path, file)
    return filename

def chipname2filename_affy(chipname):
    filename = None
    path = config.annot_data_affy
    assert os.path.exists(path),'%s does not exist'%path
    chip2file = {}
    for file in os.listdir(path):
        filename1 = os.path.join(path, file)
        chipname1,version = hash_chipname(filename1)
        if chipname1 in chip2file.keys():
            if version > chip2file[chipname1][0]:
                chip2file[chipname1] = (version,filename1)
        else:
            chip2file[chipname1] = (version,filename1)
    if chipname in chip2file.keys():
        version,filename = chip2file[chipname]
    return filename

def identify_platform_of_matrix(DATA):
    platform_list = identify_all_platforms_of_matrix(DATA)
    if not platform_list:
        return None
    out_platform = platfrom_list[0][1]
    return out_platform
    
def identify_all_platforms_of_matrix(DATA):
    """return a list of (header,platform) we can identify"""
    ids = DATA.row_names()
    chips = dict()
    for id in ids:
        x = DATA.row_names(id)
        possible_chip = identify_platform_of_annotations(x)
        if possible_chip:
            chips[possible_chip]=id
    order_platforms = arrayplatformlib.prioritize_platforms(chips.keys())
    #if chips is empty, will return an empty list
    return [(chips[platform],platform) for platform in order_platforms]

    
def identify_platform_of_annotations(annotations):
    platform2annot_cs = {}  # chip -> psid -> 1
    platform2annot_ci = {}  # chip -> psid -> 1
    paths = []
    possible_chips = []
    annotations = [i for i in annotations if len(i)>0]
    upannotations = [i.upper() for i in annotations]
    root = config.psid2platform
    assert os.path.exists(root),'the %s does not exisits'%psid2platform
    for subfolder in os.listdir(root):
        if '.DS_Store' in subfolder:
            continue
        assert os.path.isdir(os.path.join(root,subfolder))
        for platform in os.listdir(os.path.join(root,subfolder)):
            paths.append((root,subfolder,platform))
    for x in paths:
        root,subfolder,platform = x
        assert subfolder in ['case_sensitive','case_insensitive']
        f = file(os.path.join(root,subfolder,platform),'r')
        text = f.readlines()
        text = [i.strip() for i in text if len(i.strip())>0]
        f.close()
        chipname = os.path.splitext(platform)[-2] #remove the '.txt'
        
        if subfolder == 'case_insensitive':
            for psid in text:
                platform2annot_cs.setdefault(chipname, {})[psid.upper()] = 1
        else:
             for psid in text:
                platform2annot_ci.setdefault(chipname, {})[psid] = 1
                
    for chip in platform2annot_cs:
        for psid in upannotations:
            if psid not in platform2annot_cs[chip]:
                break
        else:
            possible_chips.append(chip)

    for chip in platform2annot_ci:
        for psid in annotations:
            if psid not in platform2annot_ci[chip]:
                break
        else:
            possible_chips.append(chip)  
    
    if not possible_chips:
        return None
    #combine the dict for both case_sensitive and case_insensitive into one dict
    platform2annot = platform2annot_cs.copy()
    for chip in platform2annot_ci:
        platform2annot[chip]=platform2annot_ci[chip]

    # Sort the chips by size, from smallest to largest.
    schwartz = [(len(platform2annot[chip]), chip) for chip in possible_chips]
    schwartz.sort()
    possible_chips = [x[-1] for x in schwartz]
    # Choose the smallest chip that contains all these probe sets.
    return possible_chips[0]