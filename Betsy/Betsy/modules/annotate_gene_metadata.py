#annotate_gene_metadata.py
import os
import module_utils
import subprocess
import config
import shutil
from genomicode import filelib,arrayannot
import arrayio
def run(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    outfile = get_outfile(parameters,objects,pipeline)
    slice_path = config.RENAME
    slice_BIN = module_utils.which(slice_path)
    assert slice_BIN,'cannot find the %s' %slice_path
    M = arrayio.read(single_object.identifier)
    chipname = arrayannot.identify_platform_of_matrix(M)
    if chipname:
        chipfilename = arrayannot.chipname2filename(chipname)
        if not chipfilename:
            raise AssertionError, 'I cannot annotate array platform :%s"'% chipname
        assert os.path.exists(chipfilename), 'the annotation file %s does not exist'% chipfilename
        if chipname in ['MG_U74Av2', 'HG_U133_Plus_2', 'Mu11KsubA', 'Mu11KsubB',
                        'Hu6800', 'HG_U133B', 'Mouse430_2',
                        'RG_U34A','Mouse430A_2', 'HG_U95A', 'HG_U133A', 'RAE230A',
                        'Hu35KsubC', 'Hu35KsubB', 'Hu35KsubA', 'Hu35KsubD', 'MG_U74Cv2', 'HG_U133A_2',
                        'MG_U74Bv2', 'HG_U95Av2']:
            annot_affymetrix(slice_BIN, chipfilename, single_object.identifier, outfile)
        elif chipname in ['HumanHT_12','MouseRef_8']:
            annot_illumina(slice_BIN, chipfilename, single_object.identifier, outfile)
        elif chipname in ['HumanHT_12_control','MouseRef_8_control']:
            shutil.copyfile(single_object.identifier, outfile)
        elif chipname in ['Entrez_ID_human', 'Entrez_ID_mouse', 'Entrez_symbol_human',
                          'Entrez_symbol_mouse']:
            print 'we do not annot the platform %s' % chipname
            shutil.copyfile(single_object.identifier, outfile)
        new_objects = get_newobjects(parameters, objects, pipeline)
    else:
        print 'we do not annot the platform %s' % chipname
        new_objects = objects[:]
        shutil.copyfile(single_object.identifier, outfile)
    assert module_utils.exists_nz(outfile),(
            'the output file %s for annot_file fails'% outfile)
    module_utils.write_Betsy_parameters_file(parameters,single_object,pipeline,outfile)
    return new_objects
    
def make_unique_hash(identifier,pipeline,parameters):
    return module_utils.make_unique_hash(
        identifier,pipeline,parameters)

def get_outfile(parameters,objects,pipeline):
    single_object = get_identifier(parameters,objects)
    original_file = module_utils.get_inputid(single_object.identifier)
    filename = parameters['filetype'] + '_annot_'+original_file+'.txt'
    outfile = os.path.join(os.getcwd(),filename)
    return outfile

def get_identifier(parameters,objects):
    single_object = module_utils.find_object(
        parameters,objects,parameters['filetype'],'contents,preprocess')
    assert os.path.exists(single_object.identifier),(
    'the input file %s for annot_file does not exist'%single_object.identifier)
    return single_object

def get_newobjects(parameters,objects,pipeline):
    outfile = get_outfile(parameters,objects,pipeline)
    single_object = get_identifier(parameters,objects)
    new_objects = module_utils.get_newobjects(
        outfile,parameters['filetype'],parameters,objects,single_object)
    return new_objects


def annot_affymetrix(slice_BIN,annot_file,filename,outfile):
    command1 = ['python', slice_BIN, '--remove_comments','#','--read_as_csv',
               '--clean_only',annot_file]
    f=file('annot.txt','w')
    process = subprocess.Popen(command1,shell=False,
                                stdout=f,
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
            raise ValueError(error_message)
    f.close()
    command2 = ['python',slice_BIN,'--add_row_annot',
                "annot.txt,Target Description,Entrez Gene,Gene Symbol",
                filename]
    f=file(outfile,'w')
    process = subprocess.Popen(command2,shell=False,
                                stdout=f,
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
            raise ValueError(error_message)
    f.close()

def annot_illumina(slice_BIN,annot_file,filename,outfile):
    annot_file = module_utils.gunzip(annot_file)
    f = file(annot_file,'r')
    text = f.readlines()
    start = text.index('[Probes]\n')
    end = text.index('[Controls]\n')
    f=file('annot.txt','w')
    for i in range(start+1,end):
        f.write(text[i])
    f.close()
    command2 = ['python',slice_BIN,'--add_row_annot',
                "annot.txt,Definition,Entrez_Gene_ID",
                filename]
    f=file(outfile,'w')
    process = subprocess.Popen(command2,shell=False,
                                stdout=f,
                                stderr=subprocess.PIPE)
    error_message = process.communicate()[1]
    if error_message:
            raise ValueError(error_message)
    f.close()
    
