"""
Objects:
DataObject

Functions:
get_identifier
get_inputid
make_unique_hash

write_Betsy_parameters_file
write_Betsy_report_parameters_file

which
exists_nz
is_number

convert_gene_list_platform
convert_to_same_platform

replace_matrix_header
format_convert
is_missing

plot_line_keywds
plot_line_keywd
plot_pca
find_pcaplots

download_ftp
download_dataset

gunzip
extract_from_zip
unzip_if_zip

merge_two_files
renew_parameters
high_light_path
process_group_info

"""


FMT = "%a %b %d %H:%M:%S %Y"

# Data + identifier.
# Maybe should call IdentifiedData
class DataObject:
    def __init__(self, data, identifier=""):
        self.data = data
        self.identifier = identifier
    def __repr__(self):
        x = str(self.data) + ' identifier:' + self.identifier
        return x


# TODO: don't make contents hard coded.
class AntecedentFilter:
    def __init__(self, datatype_name=None, contents=None, **attributes):
        self.datatype_name = datatype_name
        self.contents = contents
        self.attributes = attributes
    def matches_node(self, node):
        if self.datatype_name and self.datatype_name != node.datatype.name:
            return False
        if self.contents and node.attributes.get("contents") != contents:
            return False
        for (key, value) in self.attributes.iteritems():
            if node.attributes.get(key) != value:
                return False
        return True


def _find_ids_that_pass_filters(network, node_ids, filters):
    # For each of the filters, pick out one of the nodes.  Return a
    # list of nodes parallel to filters, or None if not found

    matches = []
    for f in filters:
        for id_ in node_ids:
            if id_ in matches:  # don't reuse nodes
                continue
            if f.matches_node(network.nodes[id_]):
                matches.append(id_)
                break
        else:
            return None
    assert len(matches) == len(filters)
    return matches


def find_antecedents(network, module_id, user_attributes, pool, *filters):
    # filters should be AntecedentFilter objects.  Return either a
    # DataObject (if 0 or 1 filters), or a list of DataObjects
    # parallel to filters.  Raises an exception if no antecedents
    # could be found.
    import os
    import bie3

    # Make a list of every possible combination of inputs that goes
    # into this module.
    prev_ids = []
    for id_ in network.transitions:
        if module_id in network.transitions[id_]:
            prev_ids.append(id_)
    all_input_ids = bie3._get_valid_input_combinations(
        network, module_id, prev_ids, user_attributes)

    # Filter for just the combinations in which all input nodes have
    # been run.
    filtered = []
    for input_ids in all_input_ids:
        # If not all the input nodes have been run, then ignore this
        # combination.
        x = [x for x in input_ids if x in pool]
        if len(x) != len(input_ids):
            continue
        filtered.append(x)
    all_input_ids = filtered

    # Filter based on the user criteria.
    ids = None
    if not filters:
        # If no filters given, then just return the first id found.
        assert all_input_ids
        assert all_input_ids[0]
        ids = [all_input_ids[0][0]]
    else:
        for input_ids in all_input_ids:
            ids = _find_ids_that_pass_filters(network, input_ids, filters)
            if ids is not None:
                break
    assert ids, 'cannot find node that match for %s' % \
               network.nodes[module_id].name
    for id_ in ids:
        if not pool[id_].identifier:
            continue
        assert os.path.exists(pool[id_].identifier), (
            'the input file %s for %s does not exist' %
            (pool[id_].identifier, network.nodes[module_id].name))
    objs = [pool[x] for x in ids]
    assert not filters or len(objs) == len(filters)
    if len(objs) <= 1:  # Return a single DataObject if 0 or 1 filters given.
        objs = objs[0]
    return objs

## def get_identifier(
##     network, module_id, pool, user_attributes, datatype=None, contents=None,
##     optional_key=None, optional_value=None, second_key=None, second_value=None,
##     **param):
##     # Returns a single DataObject that goes into this module.  What is
##     # this used for?
##     import os
##     import bie3

##     assert not (optional_key and not optional_value)
##     assert not (optional_value and not optional_key)
##     assert not (second_key and not second_value)
##     assert not (second_value and not second_key)

##     # Make a list of every possible combination of inputs that goes
##     # into this module.
##     prev_ids = []
##     for id_ in network.transitions:
##         if module_id in network.transitions[id_]:
##             prev_ids.append(id_)
##     all_input_ids = bie3._get_valid_input_combinations(
##         network, module_id, prev_ids, user_attributes)

##     # Filter for just the all_input_ids in which all input nodes have
##     # been run.
##     filtered = []
##     for input_ids in all_input_ids:
##         # If not all the input nodes have been run, then ignore this
##         # combination.
##         x = [x for x in input_ids if x in pool]
##         if len(x) != len(input_ids):
##             continue
##         filtered.append(x)
##     all_input_ids = filtered

##     # Filter based on the user criteria.
##     ids = []
##     for input_ids in all_input_ids:
##         for id_ in input_ids:
##             if id_ in ids:
##                 continue
##             node = network.nodes[id_]
##             if datatype and node.datatype.name != datatype:
##                 continue
##             if contents and node.attributes.get("contents") != contents:
##                 continue
##             if optional_key and \
##                    node.attributes.get(optional_key) != optional_value:
##                 continue
##             if second_key and \
##                    node.attributes.get(second_key) != second_value:
##                 continue
##             all_found = True
##             for key in param:
##                 if node.attributes.get(key) != param[key]:
##                     all_found = False
##             if not all_found:
##                 continue
##             ids.append(id_)

##     assert ids, 'cannot find node that match for %s' % \
##                network.nodes[module_id].name)
##     for id_ in ids:
##         if pool[i].identifier:
##             assert os.path.exists(pool[i].identifier), (
##                 'the input file %s for %s does not exist' %
##                 (pool[i].identifier, network.nodes[module_id].name))
##     id_ = ids[0]
##     return pool[id_]


# Not sure exactly what this does.  Looks like it takes some sort of
# file identifier, and returns some sort of cleaned up id.
def get_inputid(identifier):
    import os

    x = identifier
    x = os.path.split(x)[-1]     # get file (no path)
    x = os.path.splitext(x)[-2]  # get base name (no extension)
    x = x.split("_")[-1]         # ???
    return x
    #old_filename = os.path.split(identifier)[-1]
    #old_filename_no_ext = os.path.splitext(old_filename)[-2]
    #inputid = old_filename_no_ext.split('_')[-1]
    #return inputid


# Rename to hash something.
def make_unique_hash(identifier, pipeline, parameters, user_options):
    import os
    import hash_method

    input_file = os.path.split(identifier)[-1]
    new_parameters = parameters.copy()
    new_parameters['filesize'] = os.path.getsize(identifier)
    new_parameters['checksum'] = hash_method.get_input_checksum(identifier)
    for key in user_options:
        new_parameters[key] = user_options[key]
    hash_result = hash_method.hash_parameters(
        input_file, pipeline, **new_parameters)
    return hash_result


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
            if not rma:
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
            else:
                if (node.data.attributes['quantile_norm'] == 'yes' and
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


def merge_two_files(A_file, B_file, handle):
    """input two files and merge,write the output to handle"""
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


def which(program):
    import os

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


def write_Betsy_parameters_file(
    parameters, data_nodes, outfile, user_input, pipeline, starttime, user,
    job_name):
    import os
    import json
    import time

    f = file(os.path.join(os.getcwd(), 'Betsy_parameters.txt'), 'w')
    if isinstance(data_nodes, tuple):
        # BUG: Why is st not used?
        #st = os.stat(data_nodes[0].identifier)
        modified_time = time.strftime(FMT, time.localtime())
        text = ['Module input:', [(data_node.data.datatype.name,
                                   data_node.identifier)
                                  for data_node in data_nodes
                                 ], 'Module output parameters:', parameters,
                'Pipeline module sequence:', pipeline, 'User_input:',
                user_input, 'Start time:', starttime, 'Finish time:',
                modified_time, 'User:', user, 'Jobname:', job_name, 'Outfile:',
                outfile]
    else:
        #st = 'None'
        modified_time = 'None'
        identifier = data_nodes.identifier
        if identifier:
            #st = os.stat(identifier)
            modified_time = time.strftime(FMT, time.localtime())
        text = ['Module input:', (data_nodes.data.datatype.name,
                                  identifier), 'Module output parameters:',
                parameters, 'Pipeline module sequence:', pipeline,
                'User_input:', user_input, 'Start time:', starttime,
                'Finish time:', modified_time, 'User:', user, 'Jobname:',
                job_name, 'Outfile:', outfile]
    newtext = json.dumps(text, sort_keys=True, indent=4)
    f.write(newtext)
    f.close()


def write_Betsy_report_parameters_file(inputs, outfile, starttime, user,
                                       job_name):
    import os
    import json
    import time
    #from stat import *

    st = os.stat(outfile)
    modified_time = time.asctime(time.localtime(st[time.ST_MTIME]))
    f = file(os.path.join(os.getcwd(), 'Betsy_parameters.txt'), 'w')
    if isinstance(inputs, list):
        text = ['Module input:', [
            ('', i) for i in inputs
        ], 'Module output:', os.path.split(outfile)[-1],
                'Module output parameters:', '', 'Pipeline module sequence:',
                '', 'Start time:', starttime, 'Finish time:', modified_time,
                'User:', user, 'Jobname:', job_name]
    else:
        text = ['Module input:', ('', inputs), 'Module output:',
                os.path.split(outfile)[-1], 'Module output parameters:', '',
                'Pipeline module sequence:', '', 'Start time:', starttime,
                'Finish time:', modified_time, 'User:', user, 'Jobname:',
                job_name]
    newtext = json.dumps(text, sort_keys=True, indent=4)
    f.write(newtext)
    f.close()


def exists_nz(filename):
    """check if the filename exists and not empty"""
    import os

    if not os.path.exists(filename):  # does not exist
        return False
    if os.path.isdir(filename):  # is directory and not empty
        if os.listdir(filename):
            return True
        return False
    size = os.path.getsize(filename)  #is file and not empty
    if size > 0:
        return True
    return False


def plot_line_keywds(filename, keywords, outfile):
    import arrayio
    from genomicode import mplgraph

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
        fig = mplgraph.lineplot(*lines,
                                box_label=label,
                                legend=legend_name,
                                ylim_min=0,
                                ylabel=keyword,
                                left=0.1)
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
    import arrayio
    from genomicode import mplgraph

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
    # Bug: keywords?
    assert len(data) > 0, 'cannot find the keyword %s in the file %s' % (
        keywords, filename)
    for i in range(len(data)):
        line = [(j, data[i][j]) for j in range(len(data[i]))]
        lines.append(line)
    fig = mplgraph.lineplot(*lines,
                            box_label=label,
                            legend=legend_name,
                            ylim_min=0,
                            ylabel='Signal',
                            left=0.1)
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
    from genomicode import arrayplatformlib

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
            command = ['python', Annot_BIN, '-f', filename, '-o', 'tmp',
                       "--platform", platform]
            process = subprocess.Popen(command,
                                       shell=False,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            error_message = process.communicate()[1]
            if error_message:
                raise ValueError(error_message)
            #assert module_utils.exists_nz('tmp'), (
            #    'the platform conversion fails')
            assert exists_nz('tmp'), 'the platform conversion fails'
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
        fig = mplgraph.scatter(x1, x2,
                               xlabel='Principal Component 1',
                               ylabel='Principal Component 2',
                               color=opts,
                               legend=legend)
    else:
        fig = mplgraph.scatter(x1, x2,
                               label=labels,
                               xlabel='Principal Component 1',
                               ylabel='Principal Component 2',
                               color=opts)
    fig.savefig(result_fig)
    assert exists_nz(result_fig), 'the plot_pca.py fails'

    ##def extract_from_zip(zipName):
    ##    z = zipfile.ZipFile(zipName)
    ##    for f in z.namelist():
    ##        if f.endswith('/'):
    ##            os.makedirs(f)
    ##        else:
    ##            z.extract(f)


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


def replace_matrix_header(M, old_header, new_header):
    M = M.matrix()
    assert old_header in M._row_order
    M._row_names[new_header] = M._row_names[old_header]
    del M._row_names[old_header]
    ids = M._row_order
    ids = [x.replace(old_header, new_header) if x == old_header else x
           for x in ids]
    M._row_order = ids
    return M


def process_group_info(group_file):
    """return a dict with <sample_name:[[left_sample_list],
                                        [right_sample_list]]"""
    f = file(group_file, 'r')
    text = f.readlines()
    f.close()
    group_dict = {}
    text = [line.strip() for line in text if line.strip()]
    for line in text:
        words = line.split('\t')
        if len(words) == 3:
            if words[0] not in group_dict:
                group_dict[words[0]] = [words[2]]
            else:
                group_dict[words[0]].append(words[2])
        elif len(words) == 4:
            if words[0] not in group_dict:
                group_dict[words[0]] = [[words[2]], [words[3]]]
            else:
                group_dict[words[0]][0].append(words[2])
                group_dict[words[0]][1].append(words[3])
        else:
            raise ValueError('group file is invalid')
    return group_dict

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
