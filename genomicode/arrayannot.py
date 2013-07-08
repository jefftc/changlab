#arrayannot.py
import os
import subprocess
import tempfile
import arrayio
from genomicode import jmath, config, filelib, arrayplatformlib

def create_annot_file_affymetrix(filename):
    slice_BIN = config.slice_matrix
    command = ['python', slice_BIN, '--remove_comments', '#',
                '--read_as_csv', '--clean_only', filename]
    annot_file = None
    try:
        x,annot_file = tempfile.mkstemp(dir=".");os.close(x)
        f = file(annot_file, 'w')
        process = subprocess.Popen(command, shell=False,
                               stdout=f,
                               stderr=subprocess.PIPE)
        f.close()
        error_message = process.communicate()[1]
        if error_message:
            raise ValueError(error_message)
        matrix = [x for x in filelib.read_cols(annot_file)]
    finally:
        if annot_file and os.path.exists(annot_file):
            os.remove(annot_file)
    return matrix


def annotate_probe_affymetrix_file(probe_ids,annotation):
    headers = arrayplatformlib.affy_headers
    new_annotation = headers[annotation][0]
    platform = arrayplatformlib.identify_platform_of_annotations(probe_ids)
    filename = arrayplatformlib.chipname2filename_affy(platform)
    if not filename:
        return None
    matrix = create_annot_file_affymetrix(filename)
    header = matrix[0]
    matrix = jmath.transpose(matrix[1:])
    if new_annotation not in header:
        return None
    index = header.index('Probe Set ID')
    annot_index = header.index(new_annotation)
    result = [''] * len(probe_ids)
    for i in range(len(probe_ids)):
        gindex = matrix[index].index(probe_ids[i])
        result[i] = matrix[annot_index][gindex]
    return result


def annotate_probe_illumina_file(probe_ids, annotation):
    headers = arrayplatformlib.illu_headers
    new_annotation = headers[annotation][0]
    platform = arrayplatformlib.identify_platform_of_annotations(probe_ids)
    filename = arrayplatformlib.chipname2filename_illu(platform)
    if not filename:
        return None
    f = filelib.openfh(filename)
    text = f.readlines()
    start = text.index('[Probes]\n')
    end = text.index('[Controls]\n')
    matrix = []
    for i in range(start + 1, end):
        cols = text[i].rstrip("\r\n").split('\t')
        matrix.append(cols)
    header = matrix[0]
    matrix = jmath.transpose(matrix[1:])
    if new_annotation not in header:
        return None
    index = header.index('Probe_Id')
    annot_index = header.index(new_annotation)
    result = [''] * len(probe_ids)
    for i in range(len(probe_ids)):
        if probe_ids[i] in matrix[index]:
            gindex = matrix[index].index(probe_ids[i])
            result[i] = matrix[annot_index][gindex]
        else:
            result[i] = ' '
    return result
        

def annotate_probe_biomart(probe_ids, annotation):
    headers = arrayplatformlib.biomart_headers
    new_annotations = headers[annotation]
    platform = arrayplatformlib.identify_platform_of_annotations(probe_ids)
    in_attribute = arrayplatformlib.get_bm_attribute(platform)
    in_mart = arrayplatformlib.get_bm_organism(platform)
    R = jmath.start_R()
    jmath.R_equals_vector(probe_ids,'gene_id')
    R('library(biomaRt)')
    R('x=getOption("warn")')
    jmath.R_equals(in_attribute,'in_attribute')
    jmath.R_equals(in_attribute,'filters')
    jmath.R_equals(in_mart,'in_mart')
    R('old=useMart("ensembl",in_mart)')
    for new_annotation in new_annotations:
        try:
            jmath.R_equals(new_annotation, 'out_attribute')
            result = [''] * len(probe_ids)
            R('options(warn=-1)')
            R('homolog = getLDS(attributes=in_attribute,filters=filters,values=gene_id,mart=old,attributesL=out_attribute,martL=old)')
            break
        except:
            continue
    R('options(warn=x)')
    homolog=R['homolog']
    old_id = [str(i) for i in homolog[0]]
    human_id = [str(i) for i in homolog[1]]
    new_old_id = [old_id[i] for i in range(len(old_id)) if human_id[i]]
    new_human_id = [human_id[i] for i in range(len(human_id)) if human_id[i]]
    a = []
    for i in range(len(new_old_id)):
        gene_index = probe_ids.index(new_old_id[i])
        if gene_index in a:
            result[gene_index] = result[gene_index] + '///' + new_human_id[i] 
        else:
            result[gene_index] = new_human_id[i]
        a.append(gene_index)   
    return result


def annotate_probes_multiple(probe_ids,list_of_annotation):
    dictionary = {}
    for annotation in list_of_annotation:
        result = annotate_probes(probe_ids,annotation)
        dictionary[annotation] = result
    return dictionary


def annotate_probes(probe_ids, annotation):
    platform = arrayplatformlib.identify_platform_of_annotations(probe_ids)
    assert platform, 'we cannot identify the platform of the input probes'
    result = [''] * len(probe_ids)
    result = annotate_probe_affymetrix_file(probe_ids,annotation)
    if result is None:
        result = annotate_probe_illumina_file(probe_ids,annotation)
        if result is None:
            result = annotate_probe_biomart(probe_ids,annotation)
    return result

