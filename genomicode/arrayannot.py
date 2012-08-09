#arrayannot.py

import os
from genomicode import filelib
import Betsy_config
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
    path = Betsy_config.ANNOT_DATA_ILLU
    assert os.path.exists(path),'%s does not exist'%path
    chipname=chipname.replace('_','-')
    for file in os.listdir(path):
        if chipname in file:
            filename = os.path.join(path, file)
    return filename

def chipname2filename_affy(chipname):
    filename = None
    path = Betsy_config.ANNOT_DATA_AFFY
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



def guess_chip(filename):
    DATA = arrayio.read(filename)
    x = DATA.row_names(const.ROW_ID)
    return guess_chip_from_probesets(x)

    
def guess_chip_from_probesets(probesets):
    chip2psid1 = {}  # chip -> psid -> 1
    chip2psid2 = {}  # chip -> psid -> 1
    root = Betsy_config.PSID2PLATFORM
    paths = []
    possible_chips = []
    uprobesets = [i.upper() for i in probesets]
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
                chip2psid1.setdefault(chipname, {})[psid.upper()] = 1
        else:
             for psid in text:
                chip2psid2.setdefault(chipname, {})[psid] = 1
                
    for chip in chip2psid1:
        for psid in uprobesets:
            if psid not in chip2psid1[chip]:
                break
        else:
            possible_chips.append(chip)
    for chip in chip2psid2:
        for psid in probesets:
            if psid not in chip2psid2[chip]:
                break
        else:
            possible_chips.append(chip)  
    
    if not possible_chips:
        return None
    #combine the dict for both case_sensitive and case_insensitive into one dict
    chip2psid = chip2psid1.copy()
    for chip in chip2psid2:
        chip2psid[chip]=chip2psid2[chip]
    # Sort the chips by size, from smallest to largest.
    schwartz = [(len(chip2psid[chip]), chip) for chip in possible_chips]
    schwartz.sort()
    possible_chips = [x[-1] for x in schwartz]
    # Choose the smallest chip that contains all these probe sets.
    return possible_chips[0]
