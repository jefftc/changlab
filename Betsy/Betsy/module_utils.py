#module_utils.py
import hash_method
import arrayio
from genomicode import binreg, Matrix, jmath, matrixlib, mplgraph,arrayplatformlib
import os
import read_label_file
import json
import math
from xml.dom.minidom import parseString
import zipfile
import time
from time import strftime,localtime
from stat import *
"""contain some functions that are called by many modules"""

FMT = "%a %b %d %H:%M:%S %Y"

class Analysis:
    def __init__(self, name, parameters):
        self.name = name
        self.parameters = parameters


class DataObject:
    def __init__(self, objecttype, attributes, identifier):
        self.objecttype = objecttype
        self.attributes = attributes
        self.identifier = identifier


def get_inputid(identifier):
    old_filename = os.path.split(identifier)[-1]
    old_filename_no_ext = os.path.splitext(old_filename)[-2]
    inputid = old_filename_no_ext.split('_')[-1]
    return inputid


def make_unique_hash(identifier, pipeline, parameters):
    parameters = renew_parameters(parameters, ['status'])
    input_file = os.path.split(identifier)[-1]
    new_parameters = parameters.copy()
    new_parameters['filesize'] = os.path.getsize(identifier)
    new_parameters['checksum'] = hash_method.get_input_checksum(identifier)
    hash_result = hash_method.hash_parameters(
        input_file, pipeline, **new_parameters)
    return hash_result


def get_newobjects(outfile, out_objecttype,
                   parameters, objects, single_object):
    parameters = renew_parameters(parameters, ['status'])
    attributes = parameters.values()
    new_object = DataObject(out_objecttype, attributes, outfile)
    new_objects = objects[:]
    new_objects.remove(single_object)
    new_objects.append(new_object)
    return new_objects


def merge_two_files(A_file, B_file, handle):
    """input two files and merge,write the output to handle"""
    M_A = arrayio.read(A_file)
    M_B = arrayio.read(B_file)
    assert arrayio.tab_delimited_format.is_matrix(M_A)
    assert arrayio.tab_delimited_format.is_matrix(M_B)
    [M_A, M_B] = matrixlib.align_rows(M_A, M_B)
    assert M_A.nrow() > 0, 'there is no common genes between two files'
    X = []
    for i in range(M_A.dim()[0]):
        x = M_A._X[i] + M_B._X[i]
        X.append(x)
    row_names = M_A._row_names
    row_order = M_A._row_order
    col_names = {}
    for name in M_A._col_names:
        if name not in M_B._col_names:
            continue
        newsample_list = []
        for sample in M_B._col_names[name]:
            if sample in M_A._col_names[name]:
                newsample = sample + '_2'
            else:
                newsample = sample
            newsample_list.append(newsample)
        #x = M_A._col_names[name] + M_B._col_names[name]
        x = M_A._col_names[name] + newsample_list
        col_names[name] = x
    M_c = Matrix.InMemoryMatrix(X, row_names, col_names, row_order)
    arrayio.tab_delimited_format.write(M_c, handle)


def which(program):
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    def ext_candidates(fpath):
        yield fpath
        for ext in os.environ.get("PATHEXT", "").split(os.pathsep):
            yield fpath + ext
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            for candidate in ext_candidates(exe_file):
                if is_exe(candidate):
                    return candidate
    return None


def format_convert(X):
    data = []
    for i in range(X.dim()[1]):
        data.append(X.value(None, i))
    return data


def write_Betsy_parameters_file(parameters, single_object,
                                pipeline, outfile,starttime,user,job_name):
    st = os.stat(outfile)
    modified_time = time.asctime(time.localtime(st[ST_MTIME]))
    f = file(os.path.join(os.getcwd(), 'Betsy_parameters.txt'), 'w')
    if isinstance(single_object, list):
        text = ['Module input:', [(i.objecttype, i.identifier)
                                  for i in single_object],
                'Module output:', os.path.split(outfile)[-1],
                'Module output parameters:', parameters,
                'Pipeline module sequence:', pipeline,
                'Start time:',starttime,
                'Finish time:',modified_time,
                'User:',user,
                'Jobname:',job_name]
    else:
        text = ['Module input:', (single_object.objecttype,
                                  single_object.identifier),
                'Module output:',os.path.split(outfile)[-1],
                'Module output parameters:', parameters,
                'Pipeline module sequence:', pipeline,
                'Start time:',starttime,
                'Finish time:',modified_time,
                'User:',user,
                'Jobname:',job_name]
    newtext = json.dumps(text, sort_keys=True, indent=4)
    f.write(newtext)
    f.close()


def write_Betsy_report_parameters_file(inputs,outfile,starttime,user,job_name):
    st = os.stat(outfile)
    modified_time = time.asctime(time.localtime(st[ST_MTIME]))
    f = file(os.path.join(os.getcwd(), 'Betsy_parameters.txt'), 'w')
    if isinstance(inputs, list):
        text = ['Module input:', [('', i)
                                  for i in inputs],
                'Module output:', os.path.split(outfile)[-1],
                'Module output parameters:', '',
                'Pipeline module sequence:', '',
                'Start time:',starttime,
                'Finish time:',modified_time,
                'User:',user,
                'Jobname:',job_name]
    else:
        text = ['Module input:', ('', inputs),
                'Module output:',os.path.split(outfile)[-1],
                'Module output parameters:', '',
                'Pipeline module sequence:', '',
                'Start time:',starttime,
                'Finish time:',modified_time,
                'User:',user,
                'Jobname:',job_name]
    newtext = json.dumps(text, sort_keys=True, indent=4)
    f.write(newtext)
    f.close()


def find_object(parameters, objects, objecttype,
                attribute_key, opt_attribute=None):
    """find an object given by objecttype and attribute_key
       from a list of existing objects, if not found, return None"""
    single_object = None
    # split the compared parameter key
    attribute_keys = attribute_key.split(',')
    compare_attribute = [parameters[i] for i in attribute_keys]
    if opt_attribute:
        compare_attribute.extend(opt_attribute)
    for single_object in objects:
        attribute_match_flag = True
        if objecttype == single_object.objecttype:
            for i in compare_attribute:
                if i not in single_object.attributes:
                    attribute_match_flag = False
            if attribute_match_flag :
                return single_object
    return None


def exists_nz(filename):
    """check if the filename exists and not empty"""
    if not os.path.exists(filename): # does not exist
        return False
    if os.path.isdir(filename):  # is directory and not empty
        if os.listdir(filename):
            return True
        return False
    size = os.path.getsize(filename) #is file and not empty
    if size > 0:
       return True   
    return False



def plot_line_keywds(filename, keywords, outfile):
    M = arrayio.read(filename)
    header = M.row_names()
    label = M._col_names['_SAMPLE_NAME']
    outfiles = []
    for keyword in keywords:
        out = keyword + '.png'
        lines = []
        data = []
        legend_name = []
        for i in range(M.dim()[0]):
            if M.row_names(header[1])[i] == keyword:
                data.append(M.slice()[i])
                legend_name.append(M.row_names(header[0])[i])
        assert len(data) > 0, 'cannot find the keywords %s in the file %s' % (
            keywords, filename)
        for i in range(len(data)):
            line = [(j, data[i][j]) for j in range(len(data[i]))]
            lines.append(line)
        fig = mplgraph.lineplot(*lines, box_label=label, legend=legend_name,
                                ylim_min=0, ylabel=keyword, left=0.1)
        fig.savefig(out)
        outfiles.append(out)
    import Image
    img_w_list = []
    img_h_list = []
    imgs = []
    for i in range(len(outfiles)):
        img = Image.open(outfiles[i], 'r')
        img_w, img_h = img.size
        img_w_list.append(img_w)
        img_h_list.append(img_h)
        imgs.append(img)
    total_w = max(img_w_list) + 30
    total_h = sum(img_h_list) + 10
    background = Image.new('RGBA', (total_w, total_h), (255, 255, 255, 255))
    bg_w, bg_h = background.size
    offset_w = (bg_w - max(img_w_list)) / 2
    offset_h_list = []
    for i in range(len(img_h_list)):
        offset_h = bg_h - sum(img_h_list[i:])
        offset_h_list.append(offset_h)
    for img, offset_h in zip(imgs, offset_h_list):
        background.paste(img, (offset_w, offset_h))
    background.save(outfile)
    assert exists_nz(outfile), 'the plot_line_keywds fails'

def plot_line_keywd(filename, keyword, outfile):
    M = arrayio.read(filename)
    header = M.row_names()
    label = M._col_names['_SAMPLE_NAME']
    lines = []
    data = []
    legend_name = []
    for i in range(M.dim()[0]):
        if M.row_names(header[1])[i] == keyword:
            data.append(M.slice()[i])
            legend_name.append(keyword + '(' + M.row_names(header[0])[i] + ')')
    assert len(data) > 0, 'cannot find the keyword %s in the file %s' % (
        keywords, filename)
    for i in range(len(data)):
        line = [(j, data[i][j]) for j in range(len(data[i]))]
        lines.append(line)
    fig = mplgraph.lineplot(*lines, box_label=label, legend=legend_name,
                            ylim_min=0, ylabel='Signal', left=0.1)
    fig.savefig(outfile)
    assert exists_nz(outfile), 'the plot_line_keywd fails'

    
def renew_parameters(parameters, key_list):
    newparameters = parameters.copy()
    for key in key_list:
        if key in newparameters.keys():
            del newparameters[key]
    return newparameters


def is_number(s):
    try:
        float(s)
    except ValueError:
        return False
    return True


def download_ftp(host, path, filename):
    import ftplib
    from ftplib import FTP
    import socket
    try:
        ftp = FTP(host)
    except (socket.error, socket.gaierror), e:
        raise AssertionError('Error:cannot reach %s' % host)
    try:
        ftp.login()
    except ftplib.error_perm, e:
        ftp.quit()
        raise AssertionError('Error:cannot login anonymously')
    try:
        ftp.cwd(path)
    except ftplib.error_perm, x:
        if str(x).find('No such file') >= 0:
            raise AssertionError('cannot find the %s' % path)
    filelist = []   # to store all files
    ftp.retrlines('NLST', filelist.append)
    if filename in filelist:
        f = open(filename, 'wb')
        ftp.retrbinary('RETR ' + filename, f.write)
        f.close()
        ftp.close()
    else:
        ftp.close()
        raise AssertionError('cannot find %s in %s' % (filename, host))


def download_dataset(GSEID):
    import tarfile
    #download the tar folder from geo
    host = 'ftp.ncbi.nih.gov'
    GSE_directory = 'pub/geo/DATA/supplementary/series/' + GSEID
    download_rarfile = GSEID + '_RAW.tar'
    download_ftp(host, GSE_directory, download_rarfile)
    #untar the data folder
    GSEID_path = GSEID
    if not tarfile.is_tarfile(download_rarfile):
        raise ValueError('download file is not tar file')
    tar = tarfile.open(download_rarfile)
    tar.extractall(path=GSEID_path)
    tar.close()
    os.remove(download_rarfile)
    assert os.path.exists(GSEID_path), 'the download file %s\
                        does not exist' % GSEID_path
    assert len(os.listdir(GSEID_path)) > 0, 'the untar in \
           download_geo_dataset_GPL %s fails' % GSEID_path
    return GSEID_path


def gunzip(filename):
    import gzip
    if filename.endswith('.gz'):
        newfilename = os.path.join(
            os.getcwd(), os.path.split(os.path.splitext(filename)[0])[-1])
        #unzip the gz data
        fileObj = gzip.GzipFile(filename, 'rb')
        fileObjOut = file(newfilename, 'wb')
        while 1:
            line = fileObj.readline()
            if line == '':
                break
            fileObjOut.write(line)
        fileObj.close()
        fileObjOut.close()
        assert os.path.exists(newfilename), (
            'unzip the file %s fails' % filename)
        return newfilename
    else:
        return None

def high_light_path(network_file, pipeline, out_file):
    pipeline1 = ['start']
    pipeline1.extend(pipeline)
    f = open(network_file, 'r')
    data = f.read()
    f.close()
    dom = parseString(data)
    nodes = dom.getElementsByTagName('node')
    edges = dom.getElementsByTagName('edge')
    for analysis in pipeline1:
        for node in nodes:
            value = node.childNodes[1].getAttributeNode('value').nodeValue
            if value == analysis:
                node.childNodes[11].attributes['fill'] = '#FF33CC'
    for i in range(len(pipeline1[:-1])):
        label = pipeline1[i] + ' (pp) ' + pipeline1[i + 1]
        for edge in edges:
            if edge.childNodes[1].getAttributeNode('value').nodeValue == label:
                edge.childNodes[7].attributes['width'] = '6'
    xmlstr = dom.toxml('utf-8')
    f = open(out_file, 'w')
    f.write(xmlstr)
    f.close()


def convert_gene_list_platform(genes, platform):
    platform_list = [i.name for i in arrayplatformlib.platforms]
    assert platform in platform_list, (
        'we cannot convert to the platform %s' % platform)
    chip = arrayplatformlib.guess_chip_from_probesets(genes)
    assert chip, 'we cannot guess the platform for the input file'
    in_attribute = arrayplatformlib.get_bm_attribute(chip)
    in_mart = arrayplatformlib.get_bm_organism(chip)
    out_attribute = arrayplatformlib.get_bm_attribute(platform)
    out_mart = arrayplatformlib.get_bm_organism(platform)
    R = jmath.start_R()
    jmath.R_equals_vector(genes, 'gene_id')
    R('library(biomaRt)')
    jmath.R_equals(in_attribute, 'in_attribute')
    jmath.R_equals(in_attribute, 'filters')
    jmath.R_equals(in_mart, 'in_mart')
    R('old=useMart("ensembl",in_mart)')
    jmath.R_equals(out_attribute, 'out_attribute')
    jmath.R_equals(out_mart, 'out_mart')
    R('new=useMart("ensembl",out_mart)')
    R(str('homolog = getLDS(attributes=in_attribute,') +
      str('filters=filters,values=gene_id,mart=old,') +
      str('attributesL=out_attribute,martL=new)'))
    homolog = R['homolog']
    old_id = [str(i) for i in homolog[0]]
    human_id = [str(i) for i in homolog[1]]
    return human_id


def convert_to_same_platform(filename1, filename2, platform=None):
    M1 = arrayio.read(filename1)
    platform1 = arrayplatformlib.identify_platform_of_matrix(M1)
    M2 = arrayio.read(filename2)
    platform2 = arrayplatformlib.identify_platform_of_matrix(M2)
    if platform1 == platform2:
        return filename1, filename2
    else:
        import subprocess
        from genomicode import config
        Annot_path = config.annotate_matrix
        Annot_BIN = which(Annot_path)
        assert Annot_BIN, 'cannot find the %s' % Annot_path
        if platform1 == platform:
            filename = filename2
            newfilename1 = filename1
            newfilename2 = 'tmp'
        elif platform2 == platform:
            filename = filename1
            newfilename1 = 'tmp'
            newfilename2 = filename2
        if platform:
            command = ['python', Annot_BIN, '-f', filename,
                       '-o', 'tmp', "--platform", platform]
            process = subprocess.Popen(command, shell=False,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            error_message = process.communicate()[1]
            if error_message:
                raise ValueError(error_message)
            assert module_utils.exists_nz('tmp'), (
                'the platform conversion fails')
    return newfilename1, newfilename2


def plot_pca(filename, result_fig, opts='b', legend=None):
    from genomicode import jmath, mplgraph
    import arrayio
    R = jmath.start_R()
    jmath.R_equals(filename, 'filename')
    M = arrayio.read(filename)
    labels = M._col_names['_SAMPLE_NAME']
    data = M.slice()
    jmath.R_equals(data, 'X')
    R('NUM.COMPONENTS <- 2')
    R('S <- svd(X)')
    R('U <- S$u[,1:NUM.COMPONENTS]')
    R('D <- S$d[1:NUM.COMPONENTS]')
    # Project the data onto the first 2 components.
    R('x <- t(X) %*% U %*% diag(D)')
    x1 = R['x'][0:M.ncol()]
    x2 = R['x'][M.ncol():]
    if len(opts) > 1:
        fig = mplgraph.scatter(x1, x2, xlabel='Principal Component 1',
                               ylabel='Principal Component 2',
                               color=opts, legend=legend)
    else:
        fig = mplgraph.scatter(
            x1, x2, label=labels, xlabel='Principal Component 1',
            ylabel='Principal Component 2', color=opts)
    fig.savefig(result_fig)
    assert exists_nz(result_fig), 'the plot_pca.py fails'

    
##def extract_from_zip(zipName):
##    z = zipfile.ZipFile(zipName)
##    for f in z.namelist():
##        if f.endswith('/'):
##            os.makedirs(f)
##        else:
##            z.extract(f)

    
def extract_from_zip(zipName,outdir):
    zip = zipfile.ZipFile(zipName)
    zip.extractall(path=outdir)


def unzip_if_zip(input_name):
    if zipfile.is_zipfile(input_name):
        directory = os.path.split(input_name)[-1]
        directory = os.path.splitext(directory)[0]
        directory = os.path.join(os.getcwd(), directory)
        extract_from_zip(input_name,directory)
        for dirname in os.listdir(directory):
            if not dirname == '__MACOSX':
                directory = os.path.join(directory,dirname)
    else:
        directory = input_name
    return directory


def replace_matrix_header(M, old_header, new_header):
    M = M.matrix()
    assert old_header in M._row_order
    M._row_names[new_header] = M._row_names[old_header]
    del M._row_names[old_header]
    ids = M._row_order
    ids = [x.replace(old_header, new_header)
           if x == old_header else x for x in ids]
    M._row_order = ids
    return M
