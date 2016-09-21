"""
# Network and Nodes
get_inputid
high_light_path     Does something to network.

# Network
download_ftp
download_dataset   # move to GEO?

# Zip
gunzip
extract_from_zip
unzip_if_zip

# File operations
merge_two_files         RENAME.  Merges matrices.

# Matrix operations
is_missing              If a matrix has missing values.  Move?
format_convert          Transpose matrix?
convert_gene_list_platform
convert_to_same_platform

# Outputs
plot_line_keywds
plot_line_keywd
find_pcaplots

# Miscellaneous
is_number

find_fastq_files         Find all fastq files in a folder.
find_merged_fastq_files  Find the fastq generated by merge_reads.
find_fasta_files         Find all fasta files in a folder.
find_bam_files
read_sample_group_file
fix_sample_group_filenames
assert_sample_group_file

read_normal_cancer_file
assert_normal_cancer_samples

check_inpath
calc_max_procs_from_ram
get_physical_memory
get_num_cores
get_dirsize

get_user_option
get_config
file_exists_nz
dir_exists
splitpath
findbin
sq

read_orientation
write_orientation
read_stranded
write_stranded

"""
#FMT = "%a %b %d %H:%M:%S %Y"


# Not sure exactly what this does.  Looks like it takes some sort of
# file identifier, and returns some sort of cleaned up id.
def get_inputid(identifier):
    import os

    # identifier is in format:
    # <data type>_<file>.<ext>
    # class_label_<user_given_name>.cls
    # Pull out just the <user_given_name>.
    x = identifier
    x = os.path.split(x)[-1]     # get file (no path)
    x = os.path.splitext(x)[-2]  # get base name (no extension)
    x = x.split("_")[-1]         # ???
    return x
    #old_filename = os.path.split(identifier)[-1]
    #old_filename_no_ext = os.path.splitext(old_filename)[-2]
    #inputid = old_filename_no_ext.split('_')[-1]
    #return inputid


def high_light_path(network_file, pipeline, out_file):
    from xml.dom.minidom import parseString

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


def download_ftp(host, path, filename):
    import ftplib
    from ftplib import FTP
    import socket
    try:
        ftp = FTP(host)
    except (socket.error, socket.gaierror), e:
        raise AssertionError('Error [%s]: cannot reach %s' % (str(e), host))
    try:
        ftp.login()
    except ftplib.error_perm, e:
        ftp.quit()
        raise AssertionError('Error [%s] :cannot login anonymously' % str(e))
    try:
        ftp.cwd(path)
    except ftplib.error_perm, x:
        if str(x).find('No such file') >= 0:
            raise AssertionError('cannot find the %s' % path)
    filelist = []  # to store all files
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
    import os
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
    import os
    import gzip
    import userfile

    x = userfile._unhash_storefile(filename)
    real_name = x
    if isinstance(x, tuple):
        real_name = x[1]
    if filename.endswith('.gz') or real_name.endswith('.gz'):
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
            'unzip the file %s fails' % filename
        )
        return newfilename
    else:
        return None


def extract_from_zip(zipName, outdir):
    import zipfile

    x = zipfile.ZipFile(zipName)
    x.extractall(path=outdir)


def unzip_if_zip(input_name):
    import os
    import zipfile

    if zipfile.is_zipfile(input_name):
        directory = os.path.split(input_name)[-1]
        directory = os.path.splitext(directory)[0]
        directory = os.path.join(os.getcwd(), directory)
        extract_from_zip(input_name, directory)
        for dirname in os.listdir(directory):
            if not dirname == '__MACOSX':
                directory = os.path.join(directory, dirname)
    else:
        directory = input_name
    return directory


# Actually mergine matrices?
def merge_two_files(A_file, B_file, handle):
    """input two files and merge, write the output to handle"""
    import arrayio
    from genomicode import Matrix
    from genomicode import matrixlib

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


# Why is this here?  Should also be named has_missing_values?
def is_missing(identifier):
    import arrayio

    M = arrayio.read(identifier)
    has_missing = False
    for i in range(M.dim()[0]):
        for j in range(M.dim()[1]):
            if M._X[i][j] is None:
                has_missing = True
                break
        if has_missing:
            break
    return has_missing


def format_convert(X):
    data = []
    for i in range(X.dim()[1]):
        data.append(X.value(None, i))
    return data


def convert_gene_list_platform(genes, platform):
    from genomicode import jmath
    from genomicode import arrayplatformlib

    platform_list = [i.name for i in arrayplatformlib.platforms]
    assert platform in platform_list, (
        'we cannot convert to the platform %s' % platform
    )
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
    #old_id = [str(i) for i in homolog[0]]
    human_id = [str(i) for i in homolog[1]]
    return human_id


def convert_to_same_platform(filename1, filename2, platform=None):
    import arrayio
    import subprocess
    from genomicode import config
    from genomicode import arrayplatformlib
    from genomicode import filelib

    M1 = arrayio.read(filename1)
    platform1 = arrayplatformlib.identify_platform_of_matrix(M1)
    M2 = arrayio.read(filename2)
    platform2 = arrayplatformlib.identify_platform_of_matrix(M2)
    if platform1 == platform2:
        return filename1, filename2

    Annot_path = config.annotate_matrix
    Annot_BIN = filelib.which(Annot_path)
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
        command = [
            'python', Annot_BIN, '-f', filename, '-o', 'tmp',
            "--platform", platform]
        process = subprocess.Popen(
            command, shell=False, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
        #assert module_utils.exists_nz('tmp'), (
        #    'the platform conversion fails')
        assert filelib.exists_nz('tmp'), 'the platform conversion fails'
    return newfilename1, newfilename2


    ##def extract_from_zip(zipName):
    ##    z = zipfile.ZipFile(zipName)
    ##    for f in z.namelist():
    ##        if f.endswith('/'):
    ##            os.makedirs(f)
    ##        else:
    ##            z.extract(f)


## def replace_matrix_header(M, old_header, new_header):
##     # Actually converts from GCT to TDF format.
##     M = M.matrix()
##     assert old_header in M._row_order
##     M._row_names[new_header] = M._row_names[old_header]
##     del M._row_names[old_header]
##     ids = M._row_order
##     ids = [
##         x.replace(old_header, new_header) if x == old_header else x
##         for x in ids]
##     M._row_order = ids
##     return M


## def process_group_info(group_file):
##     """return a dict with <sample_name:[[left_sample_list],
##                                         [right_sample_list]]"""
##     f = file(group_file, 'r')
##     text = f.readlines()
##     f.close()
##     group_dict = {}
##     text = [line.strip() for line in text if line.strip()]
##     for line in text:
##         words = line.split('\t')
##         if len(words) == 3:
##             if words[0] not in group_dict:
##                 group_dict[words[0]] = [words[2]]
##             else:
##                 group_dict[words[0]].append(words[2])
##         elif len(words) == 4:
##             if words[0] not in group_dict:
##                 group_dict[words[0]] = [[words[2]], [words[3]]]
##             else:
##                 group_dict[words[0]][0].append(words[2])
##                 group_dict[words[0]][1].append(words[3])
##         else:
##             raise ValueError('group file is invalid')
##     return group_dict

##def concatenate_files(input_files,outfile):
##    with open(outfile,'w') as outfile:
##        for fname in input_files:
##            with gzip.open(fname) as infile:
##                for line in infile:
##                    outfile.write(line)
##                    
##def concatenate_multiple_line(group_dict,foldername):
##    current_dir = os.getcwd()
##    new_group_dict = {}
##    for sample_name, files in group_dict.iteritems():
##        if len(files)==1:
##            outfile = os.path.join(current_dir,sample_name+'.fastq')
##            inputfiles = [os.path.join(foldername,x) for x in files[0]]
##            concatenate_files(inputfiles,outfile)
##            new_group_dict[sample_name]=[outfile]
##        elif len(files)==2:
##            outfile_left = os.path.join(current_dir,sample_name+'_R1.fastq')
##            outfile_right = os.path.join(current_dir,sample_name+'_R2.fastq')
##            inputfiles_left = [os.path.join(foldername,x) for x in files[0]]
##            inputfiles_right = [os.path.join(foldername,x) for x in files[1]]
##    	    concatenate_files(inputfiles_left,outfile_left)
##            concatenate_files(inputfiles_right,outfile_right)
##            new_group_dict[sample_name] = [outfile_left,outfile_right]
##    return new_group_dict
##   
##           
##a=process_group_info('/home/xchen/NGS/try_RSEM/big_sample_data/sample_group.txt')
##b=concatenate_multiple_line(a,'/home/xchen/NGS/try_RSEM/big_sample_data/')
##print b


## def write_Betsy_report_parameters_file(
##     inputs, outfile, starttime, user, job_name):
##     import os
##     import json
##     import time
##     #from stat import *

##     st = os.stat(outfile)
##     modified_time = time.asctime(time.localtime(st[time.ST_MTIME]))
##     f = file(os.path.join(os.getcwd(), 'Betsy_parameters.txt'), 'w')
##     if isinstance(inputs, list):
##         text = ['Module input:', [
##             ('', i) for i in inputs
##         ], 'Module output:', os.path.split(outfile)[-1],
##                 'Module output parameters:', '', 'Pipeline module sequence:',
##                 '', 'Start time:', starttime, 'Finish time:', modified_time,
##                 'User:', user, 'Jobname:', job_name]
##     else:
##         text = ['Module input:', ('', inputs), 'Module output:',
##                 os.path.split(outfile)[-1], 'Module output parameters:', '',
##                 'Pipeline module sequence:', '', 'Start time:', starttime,
##                 'Finish time:', modified_time, 'User:', user, 'Jobname:',
##                 job_name]
##     newtext = json.dumps(text, sort_keys=True, indent=4)
##     f.write(newtext)
##     f.close()


def plot_line_keywds(filename, keywords, outfile):
    import arrayio
    from genomicode import mplgraph
    from genomicode import filelib

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
        params = {
            "box_label" : label,
            "legend" : legend_name,
            "ylim_min" : 0,
            "ylabel" : keyword,
            "left" : 0.1,
            }
        fig = mplgraph.lineplot(*lines, **params)
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
    assert filelib.exists_nz(outfile), 'the plot_line_keywds fails'


def plot_line_keywd(filename, keyword, outfile):
    import arrayio
    from genomicode import mplgraph
    from genomicode import filelib

    M = arrayio.read(filename)
    header = M.row_names()
    label = M._col_names['_SAMPLE_NAME']
    lines = []
    data = []
    legend_name = []
    for i in range(M.dim()[0]):
        if M.row_names(header[1])[i] == keyword:
            data.append(M.slice()[i])
            x = "%s (%s)" % (keyword, M.row_names(header[0])[i])
            legend_name.append(x)
    assert len(data) > 0, 'cannot find the keyword %s in the file %s' % (
        keyword, filename)
    for i in range(len(data)):
        line = [(j, data[i][j]) for j in range(len(data[i]))]
        lines.append(line)
    params = {
        "box_label" : label,
        "legend" : legend_name,
        "ylim_min" : 0,
        "ylabel" : "Signal",
        "left" : 0.1,
        }
    fig = mplgraph.lineplot(*lines, **params)
    fig.savefig(outfile)
    assert filelib.exists_nz(outfile), 'the plot_line_keywd fails'


# Why is this here?  And why is rma hard coded?
def find_pcaplots(network, pool, module_id, rma=False):
    import os

    before_pcaplot = None
    after_pcaplot = None
    for x in pool:
        node, node_id = pool[x], x
        if not node.data.datatype.name == 'PcaPlot':
            continue
        if module_id in network.transitions[node_id]:
            assert os.path.exists(node.identifier), (
                'the input file %s for %s does not exist' %
                (node.identifier, network.nodes[module_id].name))
            if (node.data.attributes['quantile_norm'] == 'no' and
                node.data.attributes['combat_norm'] == 'no' and
                node.data.attributes['shiftscale_norm'] == 'no' and
                node.data.attributes['bfrm_norm'] == 'no' and
                node.data.attributes['dwd_norm'] == 'no' and
                node.data.attributes['gene_center'] == 'no' and
                node.data.attributes['gene_normalize'] == 'no' and
                node.data.attributes['unique_genes'] == 'no' and
                node.data.attributes['platform'] == 'no' and
                node.data.attributes['duplicate_probe'] == 'no' and
                node.data.attributes['group_fc'] == 'no'):
                before_pcaplot = node
            else:
                after_pcaplot = node
        if not after_pcaplot:
            after_pcaplot = before_pcaplot
    return before_pcaplot, after_pcaplot


## def renew_parameters(parameters, key_list):
##     newparameters = parameters.copy()
##     for key in key_list:
##         if key in newparameters.keys():
##             del newparameters[key]
##     return newparameters


def is_number(s):
    try:
        float(s)
    except ValueError:
        return False
    return True


def find_fastq_files(path):
    # Return a list of the FASTQ files (full filenames) under path.
    import os
    from genomicode import filelib

    fastq_suffix = [".fq", ".fastq"]
    compress_suffix = [".gz", ".bz2", ".xz", ".zip"]
    
    all_suffixes = []
    all_suffixes.extend(fastq_suffix)
    for s1 in fastq_suffix:
        for s2 in compress_suffix:
            x = s1+s2
            all_suffixes.append(x)
    
    filenames = filelib.list_files_in_path(path)
    # Filter out the fastq files.
    i = 0
    while i < len(filenames):
        p, f = os.path.split(filenames[i])
        f_l = f.lower()
        if f.startswith("."):
            del filenames[i]
            continue
        found = False
        for s in all_suffixes:
            if f_l.endswith(s):
                found = True
                break
        if found:
            i += 1
        else:
            del filenames[i]
    return filenames


def find_merged_fastq_files(sample_group_filename, fastq_path,
                            as_dict=False, find_fasta=False):
    # Read the sample group file.  Return a list of (sample,
    # pair1.fastq, pair2.fastq).  pair2.fastq will be None for single
    # end reads.  Both files are full paths.
    # If as_dict is True, then returns a dictionary of sample ->
    # (pair1.fastq, pair2.fastq).
    import os

    assert file_exists_nz(sample_group_filename)
    assert dir_exists(fastq_path)
    sample_groups = read_sample_group_file(sample_group_filename)
    x = [x[1] for x in sample_groups]
    x = sorted({}.fromkeys(x))
    all_samples = x

    # Find the fastq files for each of the samples.
    # The merge_reads module will save fastq files in the format:
    # <fastq_path>/<sample>.fastq          # if single end
    # <fastq_path>/<sample>_<Pair>.fastq   # if paired end

    opj = os.path.join
    fastq_files = []
    for sample in all_samples:
        # Look for single or paired end fastq files.
        pair1 = pair2 = None
        se_file = opj(fastq_path, "%s.fastq" % sample)
        pe_file1 = opj(fastq_path, "%s_1.fastq" % sample)
        pe_file2 = opj(fastq_path, "%s_2.fastq" % sample)
        if find_fasta:
            se_file = opj(fastq_path, "%s.fasta" % sample)
            pe_file1 = opj(fastq_path, "%s_1.fasta" % sample)
            pe_file2 = opj(fastq_path, "%s_2.fasta" % sample)
        if find_fasta and not os.path.exists(se_file):
            se_file = opj(fastq_path, "%s.fa" % sample)
        if find_fasta and not os.path.exists(pe_file1):
            pe_file1 = opj(fastq_path, "%s_1.fa" % sample)
        if find_fasta and not os.path.exists(pe_file2):
            pe_file2 = opj(fastq_path, "%s_2.fa" % sample)
            
        if os.path.exists(se_file):
            pair1 = se_file
            assert not os.path.exists(pe_file1)
            assert not os.path.exists(pe_file2)
        elif os.path.exists(pe_file1):
            assert not os.path.exists(se_file)
            assert os.path.exists(pe_file2)
            pair1 = pe_file1
            pair2 = pe_file2
        else:
            # Cannot find the fastq file.
            similar_files = []
            for x in os.listdir(fastq_path):
                if x.find(sample) >= 0:
                    similar_files.append(x)
            x = ""
            if similar_files:
                if len(similar_files) > 5:
                    similar_files = similar_files[:5] + "..."
                if len(similar_files) > 1:
                    x = "\nSimilar files are:\n%s" % "\n".join(similar_files)
                else:
                    x = "A similar file is: %s" % similar_files[0]
            raise AssertionError, "Not found in %s: %s%s" % (
                fastq_path, sample, x)

        x = sample, pair1, pair2
        fastq_files.append(x)

    # Make sure all samples are unique.
    x1 = [x[0] for x in fastq_files]
    x2 = {}.fromkeys(x1).keys()
    assert len(x1) == len(x2), "dup sample"

    if as_dict:
        sample2fastq = {}
        for (sample, pair1, pair2) in fastq_files:
            assert sample not in sample2fastq
            sample2fastq[sample] = pair1, pair2
        fastq_files = sample2fastq
    return fastq_files


def find_fasta_files(path):
    # Return a list of the FASTA files (full filenames) under path.
    import os
    from genomicode import filelib

    fa_ext = [".fa", ".fasta"]
    compress_ext = ["", ".gz", ".xz", ".bz2"]
    EXTENSIONS = []
    for x1 in fa_ext:
        for x2 in compress_ext:
            x = "%s%s" % (x1, x2)
            EXTENSIONS.append(x)
    
    filenames = filelib.list_files_in_path(path)
    # Filter out the fastq files.
    i = 0
    while i < len(filenames):
        p, f = os.path.split(filenames[i])
        f_l = f.lower()
        if f.startswith("."):
            del filenames[i]
            continue
        for x in EXTENSIONS:
            if f_l.endswith(x):
                i += 1  # is fasta file
                break
        else:
            del filenames[i]
    return filenames


def find_sam_files(path):
    # Return a list of the .sam files (full filenames) under path.
    from genomicode import filelib

    assert dir_exists(path)
    return filelib.list_files_in_path(
        path, endswith=".sam", case_insensitive=True)


def find_bam_files(path):
    # Return a list of the .bam files (full filenames) under path.
    from genomicode import filelib

    assert dir_exists(path)
    return filelib.list_files_in_path(
        path, endswith=".bam", case_insensitive=True)


def read_sample_group_file(file_or_handle):
    # Return list of (filename, sample, pair).  pair is None, 1, or 2.
    # filename is a relative path.
    # 
    # Reads can be split across multiple files (e.g. for multiple
    # lanes), or across pairs.
    # Headers:
    # Filename  Sample  Pair
    # F1         A       1
    # F3         A       2
    # F2         A       1
    # F4         A       2
    # F5         B       1
    # F6         B       2
    #
    # - Filenames should be unique.
    # - Pair should be 1 or 2.  If single end reads, just leave blank.
    # - There can be many Filenames per Sample.  There can be many
    #   Pairs per Sample (if the reads for one pair are split).
    # - The pairs that match (1 to its 2 partner) should be next to
    #   each other in the file.
    import os
    from genomicode import filelib
    
    handle = file_or_handle
    if type(handle) is type(""):
        assert os.path.exists(file_or_handle)
        handle = filelib.openfh(handle)
        
    data = []
    for d in filelib.read_row(handle, header=1, pad_cols=""):
        assert hasattr(d, "Pair"), "Missing column: Pair"
        pair = d.Pair.strip()
        assert pair in ["", "1", "2"], "Invalid pair: %s" % d.Pair
        x = d.Filename, d.Sample, pair
        data.append(x)
        
    # Make sure filenames are unique.
    seen = {}
    for x in data:
        filename, sample, pair = x
        assert filename not in seen, "Filenames not unique: %s" % filename
        seen[filename] = 1

    # If all the Pairs are "1", then make them all blank.
    x = [x[-1] for x in data]
    x = sorted({}.fromkeys(x))
    if x == ["1"]:
        for i in range(len(data)):
            filename, sample, pair = data[i]
            data[i] = filename, sample, ""

    # For each sample, make sure there isn't a mix of paired and
    # single ended files.  It must be all single ended or all paired.
    x = [x[1] for x in data]
    all_samples = sorted({}.fromkeys(x))
    for sample in all_samples:
        x = [x[2] for x in data if x[1] == sample]
        x = sorted({}.fromkeys(x))
        if x == [""] or x == ["1"]:  # All single
            continue
        elif x == ["1", "2"]:  # All paired
            continue
        raise AssertionError, "Weird pairing: %s" % sample

    # Make sure each pair is next to each other.
    for sample in all_samples:
        pairs = [x[2] for x in data if x[1] == sample]
        # Should be all "", or a pattern of "1", "2".
        x = {}.fromkeys(pairs).keys()
        if x == [""] or x == ["1"]:  # all ""
            continue
        assert len(x) % 2 == 0, "Weird pairing: %s" % sample
        for i in range(0, len(x), 2):
            assert x[i] == "1", "Weird pairing: %s" % sample
            assert x[i+1] == "2", "Weird pairing: %s" % sample

    return data


def fix_sample_group_filenames(sample_groups, fastq_path):
    # Sometimes the file names can change, due to compression or
    # uncompression.  Try to adjust for this.
    # In sample_groups, filename is the full path.
    import os

    x = find_fastq_files(fastq_path)
    fastq_files = {}  # filename -> full path
    for x in x:
        p, f = os.path.split(x)
        fastq_files[f] = x
        
    #x = os.listdir(fastq_path)
    #fastq_files = {}.fromkeys(x)

    EXTENSIONS = [".gz", ".bz2", ".xz"]

    for i in range(len(sample_groups)):
        file_, sample, pair = sample_groups[i]

        # Try the file itself.
        to_try = [file_]
        # If the file ends with any of these extensions, then try the
        # uncompressed file.
        filestem = None
        for ext in EXTENSIONS:
            if file_.lower().endswith(ext):
                filestem = file_[:-len(ext)]
                to_try.append(filestem)
                break
        # Try the file with a different extension.
        if filestem:
            for ext in EXTENSIONS:
                x = filestem + ext
                to_try.append(x)
        # Try the file with these extensions.
        for ext in EXTENSIONS:
            x = file_ + ext
            to_try.append(x)
        fastq_filename = None
        for test_file in to_try:
            if test_file in fastq_files:
                fastq_filename = fastq_files[test_file]
                break
        if fastq_filename is None:
            # File not found.  Skip this.
            continue
        
        x = fastq_filename, sample, pair
        sample_groups[i] = x
    return sample_groups
    

def assert_sample_group_file(filename, fastq_path):
    import os
    from genomicode import filelib
    
    x = read_sample_group_file(filename)
    x = fix_sample_group_filenames(x, fastq_path)
    sample_groups = x

    # Make sure each file can be found in the fastq folder.
    # If not all found, see if it's possible the files are already merged.
    all_found = True
    filenames = [x[0] for x in sample_groups]
    for x in filenames:
        if not os.path.exists(x):
            all_found = False
            break
    if not all_found:
        possible_merged_filenames = []
        for x in sample_groups:
            filename, sample, pair = x
            ms = os.path.join(fastq_path, "%s.fastq" % sample)
            mp1 = os.path.join(fastq_path, "%s_1.fastq" % sample)
            mp2 = os.path.join(fastq_path, "%s_2.fastq" % sample)
            possible_merged_filenames += [ms, mp1, mp2]
        x = [x for x in possible_merged_filenames if os.path.exists(x)]
        if x:
            print "Fastq files may already be merged."
    filelib.assert_exists_nz_many(filenames)

    # Make sure there are no duplicate files.
    x = [x[0] for x in sample_groups]
    x = [os.path.split(x)[1] for x in x]
    files = x
    file2count = {}
    for x in files:
        file2count[x] = file2count.get(x, 0) + 1
    dups = sorted([x for x in file2count if file2count[x] > 1])
    assert not dups, "Duplicate files"
    # Make sure there are no FASTQ files without samples.
    # BROKEN. Does not account for compression, e.g. .gz.
    #all_filenames = [x[0] for x in data]
    #files = os.listdir(fastq_path)
    #for file in files:
    #    assert file in all_filenames, "Not in sample group file: %s" % file


def read_normal_cancer_file(file_or_handle):
    # Return list of (normal_sample, tumor_sample).
    import os
    from genomicode import filelib
    
    handle = file_or_handle
    if type(handle) is type(""):
        assert os.path.exists(file_or_handle)
        handle = filelib.openfh(handle)

    data = []
    for d in filelib.read_row(handle, header=1, pad_cols=""):
        assert hasattr(d, "Normal"), "Missing header: Normal"
        assert hasattr(d, "Cancer"), "Missing header: Cancer"
        ns = d.Normal
        ts = d.Cancer
        ns, ts = ns.strip(), ts.strip()
        assert ns != ts
        x = ns, ts
        data.append(x)
    return data


def assert_normal_cancer_samples(normal_cancer_data, have_samples,
                                 ignore_normal_sample=False):
    # have_samples is a list of the samples that I already have.  (Can
    # also be dictionary, set, or something that implements the "in"
    # operator.)
    from genomicode import parselib

    # Make sure files exist for all the samples.
    all_samples = []  # use a list to preserve the order of the samples.
    for (normal_sample, cancer_sample) in normal_cancer_data:
        if not ignore_normal_sample:
            if normal_sample not in all_samples:
                all_samples.append(normal_sample)
        if cancer_sample not in all_samples:
            all_samples.append(cancer_sample)
    missing = [x for x in all_samples if x not in have_samples]
    x = parselib.pretty_list(missing, max_items=5)
    s = "samples"
    if len(x) == 1:
        s = "sample"
    assert not missing, "Missing %s [%d]: %s" % (s, len(missing), x)


def check_inpath(path):
    # Rename to check_path?
    import os
    assert os.path.exists(path)
    assert os.path.isdir(path)
    return path


def calc_max_procs_from_ram(gb_per_proc, buffer=32, upper_max=None):
    # Given the number of gigabytes each processor will need,
    # calculate the maximum number of processes that should be run
    # concurrently.
    # buffer is the number of gb to leave for other processes.
    # upper_max is the absolute maximum number of processes.

    # Currently uses amount of physical memory.  Should use the
    # amount of available memory instead.
    total_bytes = get_physical_memory()
    available_bytes = get_available_memory()
    assert available_bytes <= total_bytes
    #total_gb = total_bytes/(1024*1024*1024)
    total_gb = available_bytes/(1024*1024*1024)
    max_procs = (total_gb-buffer)/gb_per_proc
    max_procs = max(max_procs, 1)
    if upper_max is not None:
        assert upper_max < 256
        max_procs = min(max_procs, upper_max)
    return max_procs


def get_available_memory():
    import psutil

    x = psutil.virtual_memory()
    return x.available


def get_physical_memory():
    # Return the amount of RAM in this machine in bytes.
    import os
    page_size = os.sysconf("SC_PAGE_SIZE")
    num_pages = os.sysconf("SC_PHYS_PAGES")
    return page_size * num_pages


def get_num_cores():
    import multiprocessing
    return multiprocessing.cpu_count()


def get_dirsize(file_or_path, followlinks=False):
    import os

    if os.path.islink(file_or_path) and not followlinks:
        return 0
    if os.path.isfile(file_or_path):
        return os.path.getsize(file_or_path)
    size = 0
    x = os.listdir(file_or_path)
    x = [os.path.join(file_or_path, x) for x in x]
    if not followlinks:
        x = [x for x in x if not os.path.islink(x)]
    x = [get_dirsize(x, followlinks=followlinks) for x in x]
    size += sum(x)
    return size

    
def get_user_option(
    user_options, name, not_empty=False, allowed_values=None, type=None,
    check_file=False):
    # not_empty means I will make sure the value is not an empty value.
    # required means the user must supply a value (even if default given).
    # allowed_values is a list of the allowed values of this option.
    # type should be a function that converts the type.
    import os
    
    assert name in user_options, "Missing option: %s" % name
    value = user_options[name]
    if not_empty:
        assert value, "Empty user option: %s" % name
    if allowed_values:
        assert value in allowed_values, "Invalid option for %s: %s" % (
            name, value)
    if value and type is not None:
        value = type(value)
    if value and check_file:
        assert os.path.exists(value), "File not found: %s" % value
    return value
    

#def make_filename(sample_filename, *args):
#    # Create a new filename using information from a sample.
#    # <PATH>/<ROOT>.<EXT>


def get_config(name, which_assert_file=False, assert_exists=False,
               quote=False):
    from genomicode import filelib
    from genomicode import config

    assert hasattr(config, name), "Not configured for genomicode: %s" % name
    x = getattr(config, name)
    if which_assert_file:
        x = filelib.which_assert(x)
    elif assert_exists:
        filelib.assert_exists(x)
    if quote:
        x = sq(x)
    return x
    


def file_exists_nz(filename):
    from genomicode import filelib
    return filelib.exists_nz(filename)


def dir_exists(path):
    import os
    if not os.path.isdir(path):
        return False
    if not os.path.exists(path):
        return False
    return True


def root2filename(filenames):
    # filenames is a list of <filename>s in format:
    # <directory>/<root><ext>
    # Return a dictionary of <root> -> <filename>
    root2filename = {}
    for filename in filenames:
        path, root, ext = splitpath(filename)
        root2filename[root] = filename
    return root2filename


def splitpath(path):
    # Return tuple of <directory>, <root>, <ext>.
    # <directory>/<root><ext>
    # <ext> includes the dot.  If there is no extension (no dot), then
    # <ext> will be an empty string.
    # 
    # Examples:
    # path           directory  root  ext
    # /etc/rc.conf   /etc       rc    .conf
    # /etc/rc.conf   /etc       rc
    # /usr/etc       /usr       etc
    # /usr/etc/      /usr/etc
    # /              /
    import os
    
    dir_, f = os.path.split(path)
    root, ext = os.path.splitext(f)
    return dir_, root, ext


def findbin(name, quote=False):
    # name is the name in the genomicode config file.
    # DEPRECATE THIS FUNCTION.  Just use get_config
    return get_config(name, which_assert_file=True, quote=quote)


def sq(name):
    # quote for a shell command.
    from genomicode import parallel
    return parallel.quote(name)
    

class Orientation:
    def __init__(self, orientation, reads_ns, reads_fr, reads_rf, reads_ff):
        assert orientation in ["single", "paired_fr", "paired_rf", "paired_ff"]
        self.orientation = orientation
        self.reads_ns = reads_ns
        self.reads_fr = reads_fr
        self.reads_rf = reads_rf
        self.reads_ff = reads_ff


def read_orientation(filename):
    import json
    from genomicode import filelib
    filelib.assert_exists_nz(filename)
    text = open(filename).read()
    x = json.loads(text)
    return Orientation(**x)


def write_orientation(orient, filename):
    import json
    handle = open(filename, 'w')
    members = ["orientation", "reads_ns", "reads_fr", "reads_rf", "reads_ff"]
    x = {}
    for m in members:
        x[m] = getattr(orient, m)
    json.dump(x, handle, indent=2)


class Stranded:
    def __init__(self, single_or_paired, stranded, frac_failed, frac_first,
                 frac_second):
        assert single_or_paired in ["single", "paired"]
        assert stranded in ["unstranded", "firststrand", "secondstrand"]
        self.single_or_paired = single_or_paired
        self.stranded = stranded
        self.frac_failed = frac_failed
        self.frac_first = frac_first
        self.frac_second = frac_second


def read_stranded(filename):
    import json
    from genomicode import filelib
    filelib.assert_exists_nz(filename)
    text = open(filename).read()
    x = json.loads(text)
    return Stranded(**x)


def write_stranded(stranded, filename):
    import json
    handle = open(filename, 'w')
    members = [
        "single_or_paired", "stranded",
        "frac_failed", "frac_first", "frac_second"]
    x = {}
    for m in members:
        x[m] = getattr(stranded, m)
    json.dump(x, handle, indent=2)
