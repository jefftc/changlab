# convert_gene_ids
# 
# _convert_gene_ids_biomart
# _convert_gene_ids_local
# _convert_entrez_symbol_to_entrez
# _find_entrez_gene
#
# _clean_id
# _clean_genes_for_biomart
# _remove_dups
# _remove_refseq_version
# _start_R

GLOBAL_R = None
def _start_R():
    global GLOBAL_R
    from genomicode import jmath

    if GLOBAL_R is None:
        R = jmath.start_R()
        R('library(biomaRt)')
        GLOBAL_R = R
    return GLOBAL_R


def convert_gene_ids(
    gene_ids, in_platform, out_platform, in_delim, out_delim,
    keep_dups, keep_emptys, no_na):
    # Return a list of the output IDs, parallel to gene_ids.  If a
    # gene ID to multiple output IDs, the output IDs will be separated
    # by out_delim.  If it is missing, then the output will be an
    # empty string.
    # in_platform and out_platform are the names of the platforms.

    # Make a cleaned up version of the gene_ids to convert.
    x = []
    for gene_id in gene_ids:
        x.extend(_clean_id(gene_id, in_delim, in_platform))
    # No duplicates.
    x = {}.fromkeys(x).keys()
    gene_ids_c = x
    
    in2out = None
    if in_platform == "Entrez_Symbol_human" and \
       out_platform in ("Entrez_Symbol_human", "Entrez_ID_human"):
        in2out = _convert_entrez_symbol_to_entrez(
            gene_ids_c, out_platform, "9606")
    if in2out is None:
        in2out = _convert_gene_ids_local(in_platform, out_platform)
    if in2out is None:
        in2out = _convert_gene_ids_biomart(
            gene_ids_c, in_platform, out_platform, no_na)
    assert in2out, "I could not convert %s to %s" % (in_platform, out_platform)

    # Make a parallel list of the output IDs.
    output_ids = []
    for gene_id in gene_ids:
        in_ids = _clean_id(gene_id, in_delim, in_platform)
        out_ids = []
        for x in in_ids:
            y = in2out.get(x, [""])
            out_ids.extend(y)
        if not keep_emptys:
            out_ids = [x for x in out_ids if x]
        if not keep_dups:
            out_ids = _remove_dups(out_ids)
        x = out_delim.join(out_ids)
        output_ids.append(x)
    return output_ids


def _clean_genes_for_biomart(gene_ids):
    gene_ids = [i.replace("'",'') for i in gene_ids] #modify "2'-PDE"
    return gene_ids


def _convert_gene_ids_biomart(gene_ids, in_platform, out_platform, no_na):
    # Maximum number of genes to request at a time.
    MAX_GENES = 25000

    in2out = {}
    while gene_ids:
        batch = gene_ids[:MAX_GENES]
        gene_ids = gene_ids[MAX_GENES:]
        
        x = _convert_gene_ids_biomart_h(
            batch, in_platform, out_platform, no_na)
        in2out.update(x)
    return in2out


def _convert_gene_ids_biomart_h(gene_ids, in_platform, out_platform, no_na):
    # Return a dictionary of gene_id -> list of converted_ids, or None
    # if these platforms cannot be converted.
    from genomicode import jmath
    from genomicode import arrayplatformlib

    if not gene_ids:
        return {}
    R_fn, R_var = jmath.R_fn, jmath.R_var

    # An attribute is the biomart name for the platform.
    in_attribute = arrayplatformlib.get_bm_attribute(in_platform)
    out_attribute = arrayplatformlib.get_bm_attribute(out_platform)
    #assert in_attribute, "Bad platform: %s" % in_platform
    #assert out_attribute, "Bad platform: %s" % out_platform
    if not in_attribute or not out_attribute:
        return None
    in_mart = arrayplatformlib.get_bm_organism(in_platform)
    out_mart = arrayplatformlib.get_bm_organism(out_platform)
    assert in_mart, "No bm organism for platform: %s" % in_platform

    R = _start_R()
    gene_ids = _clean_genes_for_biomart(gene_ids)
    jmath.R_equals_vector(gene_ids, 'gene.ids')
    # Select the BioMart dataset to use.
    #mart = "ensembl"
    mart = "ENSEMBL_MART_ENSEMBL"  # Changed 151120.
    host = "www.ensembl.org"
    R_fn("useMart", mart, in_mart, host=host, RETVAL="in.dataset")
    R_fn("useMart", mart, out_mart, host=host, RETVAL="out.dataset")

    # Link two data sets and retrieve information from the linked datasets.
    R_fn(
        "getLDS", attributes=in_attribute, filters=in_attribute,
        values=R_var("gene.ids"), mart=R_var("in.dataset"),
        attributesL=out_attribute, martL=R_var("out.dataset"),
        RETVAL="homolog")

    homolog = R["homolog"]
    # homolog is DataFrame with two parallel rows:
    # <in_ids>
    # <out_ids>
    assert len(homolog) == 2, \
           "BioMart returned no results mapping from %s:%s to %s:%s." % (
        in_mart, in_attribute, out_mart, out_attribute)

    in_ids = [str(x) for x in homolog[0]]
    out_ids = [str(x) for x in homolog[1]]

    # Sometimes BioMart will generate "NA" if something is missing.
    if no_na:
        for i in range(len(out_ids)):
            if out_ids[i].upper() == "NA":
                out_ids[i] = ""
    
    in2out = {}
    for x, y in zip(in_ids, out_ids):
        if not y.strip():
            continue
        val = in2out.get(x, [])
        val.append(y)
        in2out[x] = sorted(val)

    return in2out


def _convert_gene_ids_local(in_platform, out_platform):
    # Return a dictionary of gene_id -> list of converted_ids, or None
    # if these platforms cannot be converted.
    import os
    from genomicode import config
    from genomicode import filelib

    assert os.path.exists(config.convert_platform)
    x = "%s___%s.txt" % (in_platform, out_platform)
    filename = os.path.join(config.convert_platform, x)
    if not os.path.exists(filename):
        return None

    in2out = {}
    for cols in filelib.read_cols(filename):
        # <in_id>  <out_id1> ... <out_idn>
        assert len(cols) >= 2
        in_id = cols[0]
        out_ids = cols[1:]
        in2out[in_id] = out_ids
    return in2out
        

FOUND_ID2ENTREZ = {}

def _convert_entrez_symbol_to_entrez(gene_ids, out_platform, db_tax_id):
    # Return a dictionary of gene_id -> list of converted_ids.
    # Need db_tax_id to limit search space.
    global FOUND_ID2ENTREZ

    FOUND = FOUND_ID2ENTREZ

    # Not sure if other organisms are implemented.
    assert db_tax_id in ["9606", "10090"]
    #db_tax_id = "9606"
    #db_tax_id = "10090"

    in2out = {}

    for in_id in gene_ids:
        gene_id = gene_symbol = None

        # Try to find the gene ID from the symbol.
        # First, look to see if it's already been found.
        if in_id in FOUND:
            gene_id, gene_symbol = FOUND[in_id]
        if not gene_id:
            x = _find_entrez_gene(in_id, db_tax_id)
            if x:
                gene_id, gene_symbol = x

        if out_platform == "Entrez_ID_human" and gene_id:
            in2out[in_id] = [str(gene_id)]
        elif out_platform == "Entrez_Symbol_human" and gene_symbol:
            in2out[in_id] = [gene_symbol]
    return in2out


FIND_GENE_ERROR = None
def _find_entrez_gene(gene_symbol, tax_id):
    # Return tuple of (gene_id, gene_symbol) or None.
    global FIND_GENE_ERROR
    from genomicode import genefinder

    FIND_GENE_ERROR = None

    # Try to find the gene.
    try:
        x = genefinder.find_gene(gene_symbol, tax_id=tax_id)
    except AssertionError, x:
        FIND_GENE_ERROR = str(x)
        return None
    if x:
        gene_id = x[0]
        gene_symbol = x[1]
        return gene_id, gene_symbol

    # Could not find gene.  Look for a discontinued gene.
    try:
        x = genefinder.find_discontinued(gene_symbol, tax_id=tax_id)
        if x:
            x = genefinder.find_gene(x, tax_id=tax_id)
    except AssertionError, x:
        FIND_GENE_ERROR = str(x)
        return None
    if x:
        gene_id = x[0]
        gene_symbol = x[1]
        return gene_id, gene_symbol

    return None
    

def _remove_dups(ids):
    ids_c = []
    for x in ids:
        if x not in ids_c:
            ids_c.append(x)
    return ids_c

def _remove_refseq_version(refseq_id):
    i = refseq_id.find(".")
    if i < 0:
        return refseq_id
    return refseq_id[:i]


def _clean_id(gene_id, delimiter, in_platform):
    # Return a list of this ID cleaned up.  gene_id might be split
    # into multiple IDs by the delimiter.
    x = [gene_id]
    if delimiter:
        x = gene_id.split(delimiter)
    # Clean up whitespace.
    x = [x.strip() for x in x]
    # No empty IDs.
    x = [x for x in x if x]
    # Ignore "---".
    x = [x for x in x if x != "---"]
    # Hack: Remove version numbers from RefSeq IDs.
    if in_platform.lower().find("refseq") >= 0:
        x = [_remove_refseq_version(x) for x in x]
    # No duplicates.
    x = {}.fromkeys(x).keys()
    return x



## """

## Functions:
## create_annot_file_affymetrix
## annotate_probe_affymetrix_file
## annotate_probe_illumina_file
## annotate_probe_biomart

## annotate_probes_multiple
## annotate_probes

## read_mapping_file
## map_probes

## convert_probe_ids
## _convert_probe_ids_local

## """
## import os
## import subprocess
## import tempfile
## import arrayio
## from genomicode import jmath, config, filelib, arrayplatformlib
## import re


## def create_annot_file_affymetrix(filename):
##     slice_BIN = config.slice_matrix
##     command = ['python', slice_BIN, '--remove_comments', '#',
##                 '--read_as_csv', '--clean_only', filename]
##     annot_file = None
##     try:
##         x,annot_file = tempfile.mkstemp(dir=".");os.close(x)
##         f = file(annot_file, 'w')
##         process = subprocess.Popen(command, shell=False,
##                                stdout=f,
##                                stderr=subprocess.PIPE)
##         f.close()
##         error_message = process.communicate()[1]
##         if error_message:
##             raise ValueError(error_message)
##         matrix = [x for x in filelib.read_cols(annot_file)]
##     finally:
##         if annot_file and os.path.exists(annot_file):
##             os.remove(annot_file)
##     return matrix


## def annotate_probe_affymetrix_file(probe_ids,annotation):
##     headers = arrayplatformlib.affy_headers
##     new_annotation = headers[annotation][0]
##     platform = arrayplatformlib.identify_platform_of_annotations(probe_ids)
##     filename = arrayplatformlib.chipname2filename_affy(platform)
##     if not filename:
##         return None
##     matrix = create_annot_file_affymetrix(filename)
##     header = matrix[0]
##     matrix = jmath.transpose(matrix[1:])
##     if new_annotation not in header:
##         return None
##     index = header.index('Probe Set ID')
##     annot_index = header.index(new_annotation)
##     result = [''] * len(probe_ids)
##     for i in range(len(probe_ids)):
##         gindex = matrix[index].index(probe_ids[i])
##         result[i] = matrix[annot_index][gindex]
##     return result


## def annotate_probe_illumina_file(probe_ids, annotation):
##     headers = arrayplatformlib.illu_headers
##     new_annotation = headers[annotation][0]
##     platform = arrayplatformlib.identify_platform_of_annotations(probe_ids)
##     filename = arrayplatformlib.chipname2filename_illu(platform)
##     if not filename:
##         return None
##     f = filelib.openfh(filename)
##     text = f.readlines()
##     start = text.index('[Probes]\n')
##     end = text.index('[Controls]\n')
##     matrix = []
##     for i in range(start + 1, end):
##         cols = text[i].rstrip("\r\n").split('\t')
##         matrix.append(cols)
##     header = matrix[0]
##     matrix = jmath.transpose(matrix[1:])
##     if new_annotation not in header:
##         return None
##     index = header.index('Probe_Id')
##     annot_index = header.index(new_annotation)
##     result = [''] * len(probe_ids)
##     for i in range(len(probe_ids)):
##         if probe_ids[i] in matrix[index]:
##             gindex = matrix[index].index(probe_ids[i])
##             result[i] = matrix[annot_index][gindex]
##         else:
##             result[i] = ' '
##     return result
        

## def annotate_probe_biomart(probe_ids, annotation):
##     headers = arrayplatformlib.biomart_headers
##     new_annotations = headers[annotation]
##     platform = arrayplatformlib.identify_platform_of_annotations(probe_ids)
##     in_attribute = arrayplatformlib.get_bm_attribute(platform)
##     in_mart = arrayplatformlib.get_bm_organism(platform)
##     assert in_attribute is not None, "No BioMart attribute: %s" % platform
##     assert in_mart is not None, "No BioMart mart: %s" % platform
    
##     R = jmath.start_R()
##     jmath.R_equals_vector(probe_ids,'gene_id')
##     R('library(biomaRt)')
##     R('x=getOption("warn")')
##     jmath.R_equals(in_attribute,'in_attribute')
##     jmath.R_equals(in_attribute,'filters')
##     jmath.R_equals(in_mart,'in_mart')
##     R('old=useMart("ensembl",in_mart)')
##     for new_annotation in new_annotations:
##         try:
##             jmath.R_equals(new_annotation, 'out_attribute')
##             result = [''] * len(probe_ids)
##             R('options(warn=-1)')
##             R('homolog = getLDS(attributes=in_attribute,filters=filters,values=gene_id,mart=old,attributesL=out_attribute,martL=old)')
##             break
##         except:
##             continue
##     R('options(warn=x)')
##     homolog=R['homolog']
##     old_id = [str(i) for i in homolog[0]]
##     human_id = [str(i) for i in homolog[1]]
##     new_old_id = [old_id[i] for i in range(len(old_id)) if human_id[i]]
##     new_human_id = [human_id[i] for i in range(len(human_id)) if human_id[i]]
##     a = []
##     for i in range(len(new_old_id)):
##         gene_index = probe_ids.index(new_old_id[i])
##         if gene_index in a:
##             result[gene_index] = result[gene_index] + '///' + new_human_id[i] 
##         else:
##             result[gene_index] = new_human_id[i]
##         a.append(gene_index)   
##     return result


## def annotate_probes_multiple(probe_ids, list_of_annotation):
##     dictionary = {}
##     for annotation in list_of_annotation:
##         result = annotate_probes(probe_ids, annotation)
##         dictionary[annotation] = result
##     return dictionary


## def annotate_probes(probe_ids, annotation):
##     platform = arrayplatformlib.identify_platform_of_annotations(probe_ids)
##     assert platform, 'we cannot identify the platform of the input probes'
##     result = [''] * len(probe_ids)
##     result = annotate_probe_affymetrix_file(probe_ids, annotation)
##     if result is None:
##         result = annotate_probe_illumina_file(probe_ids, annotation)
##     if result is None:
##         result = annotate_probe_biomart(probe_ids, annotation)
##     return result

## def read_mapping_file(mapping_file):
##     """read the platform file and return a mapping 
##     dictionary <probe_id: a list of new_probe_id>"""
##     probe_map = {}
##     f = file(mapping_file,'r')
##     text = f.read()
##     line1 = re.split('\r|\n', text)
##     f.close()
##     for line in line1[1:]:
##         words = line.split('\t')
##         if len(words) == 1:
##             words.append('')
##         if words[0] not in probe_map:
##             probe_map[words[0]] = [words[1]]
##         else:
##             probe_map[words[0]] = probe_map[words[0]].append(words[1])
##     return probe_map
    
## def map_probes(probe_ids, mapping_file, in_delim='///', out_delim='///'):
##     """given a list of probe_ids and the platform_file for 
##     mapping, return a list of probe_ids in new platform"""
##     probe_map = read_mapping_file(mapping_file)
##     new_probes = []
##     for probe_id in probe_ids:
##         multiple_in_ids = probe_id.split(in_delim)
##         multiple_out_ids = [probe_map.get(x) for x in multiple_in_ids]
##         multiple_out_ids = [x for x in multiple_out_ids if x]
##         multiple_out_ids = sum(multiple_out_ids,[])
##         multiple_out_ids = sorted({}.fromkeys(multiple_out_ids))
##         newid = out_delim.join(multiple_out_ids)
##         new_probes.append(newid)
##     return new_probes


## def convert_probe_ids(probe_ids, platform_name):
##     new_probes = _convert_probe_ids_local(probe_ids, platform_name)
##     if not new_probes:
##        new_probes = annotate_probe_biomart(probe_ids, platform_name)
##     return new_probes


## def _convert_probe_ids_local(probe_ids, platform_name):
##     platform_path = config.convert_platform
##     assert os.path.exists(platform_path)
##     old_platform = arrayplatformlib.identify_platform_of_annotations(probe_ids)
##     filename = old_platform + '___' + platform_name + '.txt'
##     if not os.path.exists(os.path.join(platform_path,filename)):
##         return None
##     new_probes = map_probes(probe_ids, os.path.join(platform_path,filename))
##     return new_probes


