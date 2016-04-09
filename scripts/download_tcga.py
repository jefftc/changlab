#!/usr/bin/env python

import urllib2
import re
import argparse
import os
from genomicode import parselib, filelib, arrayannot
import arrayio
from genomicode import genefinder, timer, Matrix, matrixlib

# retrieve_all_dates
# retrieve_diseases

# Not quite a good idea.
datatype_match = {
    'RSEM_genes' : (
        'Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__' +
        'RSEM_genes_normalized__data.Level_3'),
    'RSEM_exons' : 'exon_expression__data.Level_3',
    'humanmethylation450' : (
        'Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__' +
        'within_bioassay_data_set_function__data.Level_3'),
    'mirnaseq' : (
            'Merge_mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__' +
            'miR_gene_expression__data.Level_3'),
    'clinical' : 'Merge_Clinical.Level_1',
    'rppa' : '.RPPA_AnnotateWithGene.Level_3',
    'RSEM_isoforms' : (
            'Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__' +
            'RSEM_isoforms_normalized__data.Level_3'),
    'cnv_gistic2' : 'CopyNumber_Gistic2.Level_4',
    
    'agilent' : [
        'Merge_transcriptome__agilentg4502a_07_1__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.Level_3',
        'Merge_transcriptome__agilentg4502a_07_2__unc_edu__Level_3__unc_lowess_normalization_gene_level__data.Level_3',
        ],
    'affymetrix' : (
            'Merge_transcriptome__ht_hg_u133a__broad_mit_edu__Level_3__' +
            'gene_rma__data.Level_3'),
    }

datatype2resource = {
    'RSEM_genes' : 'stddata',
    'RSEM_exons' : 'stddata',
    'humanmethylation450' : 'stddata',
    'mirnaseq' : 'stddata',
    'clinical' : 'stddata',
    'rppa' : 'stddata',
    'RSEM_isoforms' : 'stddata',
    'cnv_gistic2' : 'analyses',
    'agilent' : 'stddata',
    'affymetrix' : 'stddata',
    }

resources = ['stddata', 'analyses']

URL2HTML = {}
def read_url(url):
    global URL2HTML
    if url not in URL2HTML:
        #print "Reading %s" % url; sys.stdout.flush()
        timer.wait(2, 'tcga')
        response = urllib2.urlopen(url)
        html = response.read()
        URL2HTML[url] = html
    return URL2HTML[url]

def retieve_dates_from_resource(resource):
    #import pprint  ###
    # Return a dictionary of the date -> url with given resource.
    #url points to an HTML page for the runs for a specific date.  It has a table of the
    # diseases and HREFs to the data for the diseases at that date.
    #download_url = None
    url = None
    if resource == 'stddata':
        dashboard_url = (
        'https://confluence.broadinstitute.org/display/GDAC/Dashboard-Stddata')
        url = 'http://gdac.broadinstitute.org/runs/info/stddata__runs_list.html'
    elif resource == 'analyses':
        dashboard_url =(
            'https://confluence.broadinstitute.org/display/GDAC/Dashboard-Analyses')
        url = 'http://gdac.broadinstitute.org/runs/info/analyses__runs_list.html'
    else:
        raise ValueError('we do not recognize the resource %s' %resource)
    # The latest data is listed on the "Dashboard" page.  The old data
    # is shown on a different one.  Download them separately.
    all_dates = {}
    urls = read_and_extract_urls(dashboard_url)
    diseases, date = get_disease_lastest(urls,resource)
    date = date.replace("_", "")
    all_dates[date] = dashboard_url

    # Read the old data.
    html = read_url(url)
    link_tags = parselib.get_tags_and_contents(html,'a')
    for i in link_tags:
        url = parselib.get_href(i)
        date = re.search(r'[0-9]{4}_[0-9]{2}_[0-9]{2}', url)
        if not date:
            continue
        date = date.group(0)
        date = date.replace("_", "")
        all_dates[date] = url
    assert all_dates, "No dates found"
    return all_dates


def retrieve_all_dates():
    # Return a dictionary of the resource-> date and date is a dict which date -> url.
    # url points to an HTML page for the runs for a specific date.  It has a table of the
    # diseases and HREFs to the data for the diseases at that date.
    resources = ['stddata','analyses']
    all_dates = {}  # resource->date
    all_dates_list = []
    for resource in resources:
        all_date_from_resource = retieve_dates_from_resource(resource)
        all_dates[resource] = all_date_from_resource
        all_dates_list.extend(all_date_from_resource.keys())
    all_dates_list = list(set(all_dates_list))
    
    all_dates_list = sorted(all_dates_list)
    return all_dates, all_dates_list


def extract_all_hrefs(html):
    x = parselib.get_tags_and_contents(html,'a')
    x = [parselib.get_href(x) for x in x]
    x = [x for x in x if x]
    x = [x for x in x if not x.startswith("#")]
    return x

def retrieve_disease_resource(date,resource):
    all_dates_from_resource = retieve_dates_from_resource(resource)
    if date not in all_dates_from_resource:
        return None
    url = all_dates_from_resource[date]
    # URL:
    # http://gdac.broadinstitute.org/runs/stddata__2014_07_15/data/ACC/20140715
    pattern = re.compile(
        r'http://%s/runs/%s__[0-9_]{10}/data/([A-Z]+)/([0-9]{8})' % (
            "gdac.broadinstitute.org",resource), re.IGNORECASE)
    diseases = []
    html = read_url(url)
    for href in extract_all_hrefs(html):
        m = pattern.match(href)
        if not m:
            continue
        x = m.group(1)
        diseases.append(x)
    return diseases

def retrieve_diseases(date):
    all_diseases = []
    for resource in resources:
        disease = retrieve_disease_resource(date,resource)
        if disease:
            all_diseases.extend(disease)
    all_diseases = list(set(all_diseases))
    assert all_diseases, "could not find diseases"
    return all_diseases


def download_file(disease, date, datatype):
    assert len(date) == 8
    long_date = "%s_%s_%s" % (date[:4], date[4:6], date[6:8])
    resource = datatype2resource[datatype]
    # URL:
    # http://gdac.broadinstitute.org/runs/stddata__2014_07_15/data/ACC/20140715
    link = "http://%s/runs/%s__%s/data/%s/%s/" % (
        "gdac.broadinstitute.org", resource, long_date, disease, date)
    newlinks = get_all_datas_on_page(link)
    if isinstance(datatype_match[datatype],str):
        match_items = [datatype_match[datatype]]
    elif isinstance(datatype_match[datatype],list):
        match_items = datatype_match[datatype]
    else:
        raise ValueError('cannot recognize datatype in datatype_match dict')
    resultlinks = []
    for newlink in newlinks:
        for match_item in match_items:
            if match_item in newlink and disease+'-FFPE' not in newlink:
                data = read_url(link+newlink)
                with open(newlink, "wb") as code:
                     code.write(data)
                print 'finished download %s' % newlink
                resultlinks.append(newlink)
    if resultlinks:
        return resultlinks    
    assert ValueError('download fails')

def get_data_type_resource(disease,date,resource):
    long_date = "%s_%s_%s" % (date[:4], date[4:6], date[6:8])
    # URL:
    # http://gdac.broadinstitute.org/runs/stddata__2014_07_15/data/ACC/20140715
    link = "http://%s/runs/%s__%s/data/%s/%s/" % (
        "gdac.broadinstitute.org",resource, long_date, disease, date)
    newlinks = get_all_datas_on_page(link)
    result = []
    for newlink in newlinks:
        for datatype in datatype_match:
            if datatype2resource[datatype] != resource:
                continue
            if datatype_match[datatype] not in newlink:
                continue
            # Brittle
            if disease+'-FFPE' in newlink:
                continue
            assert datatype not in result, "dup datatype"
            result.append(datatype)
    return result
    
def get_data_type(disease, date):
    assert len(date) == 8
    #long_date = "%s_%s_%s" % (date[:4], date[4:6], date[6:8])
    result = []
    for resource in resources:
       result.extend(get_data_type_resource(disease,date,resource))
    return result
                
            
def get_disease_lastest(urls,resource):
    # Return list of diseases, latest_date
    diseases = []
    pattern = re.compile(
        r'http://gdac.broadinstitute.org/runs/%s__[0-9_]*/[A-Z]*.html'%resource)
    lastest = None
    for url in urls:
        disease_link = re.findall(pattern, url)
        for link in disease_link:
            x = re.search(r'[A-Z]+', link)
            if not x:
                continue
            disease = x.group(0)
            assert disease not in diseases
            diseases.append(disease)

            x = re.search(r'[0-9]{4}_[0-9]{2}_[0-9]{2}', link)
            assert x, "date not found"
            date = x.group(0)
            assert not lastest or lastest == date
            lastest = date
    return diseases, lastest 


def get_all_datas_on_page(page):
    html = read_url(page)
    links = []
    link_tags = parselib.get_tags_and_contents(html,'a')
    for i in link_tags:
        url = parselib.get_href(i)
        # This is brittle.
        if url.startswith("gdac.broadinstitute.org_") and url.endswith(".gz"):
            links.append(url)
    return links


def read_and_extract_urls(page):
    # Read a webpage a return a list of the URLs on that page.
    html = read_url(page)
    x = parselib.get_tags_and_contents(html,'a')
    x = [parselib.get_href(x) for x in x]
    x = [x for x in x if x]
    return x

def merge_files(input_list, outfile):
    """input two files and merge,write to the outfile"""
    assert len(input_list) == 2
    A_file = input_list[0]
    B_file = input_list[1]
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
    handle = file(outfile,'w')
    arrayio.tab_delimited_format.write(M_c, handle)
    handle.close()
    
def extract_and_merge_files(gzfile_list, resource):
    result = []
    for gzfile in gzfile_list:
        if not gzfile.endswith('tar.gz'):
            return gzfile
        import tarfile
        tfile = tarfile.open(gzfile, 'r:gz')
        gzname = os.path.split(gzfile)[-1]
        newdir = os.path.join(os.getcwd(),gzname[:-7])
        tfile.extractall(newdir)
        folders = os.listdir(newdir)
        folder = None
        for folder in folders:
            if folder == '._.DS_Store':
                continue
        directory = os.path.join(newdir,folder)
        assert os.path.exists(directory)
        files = os.listdir(directory)
        if resource == 'stddata':
            for filename in files:
                if filename.endswith('txt') and filename != 'MANIFEST.txt':
                    result.append(os.path.join(directory,filename))
        elif resource == 'analyses':
            for filename in files:
                if filename =='all_data_by_genes.txt':
                    result.append(os.path.join(directory,filename))
        else:
            raise ValueError('not recoginzed resource %s'%resource)
    if len(result) == 1:
        return result[0]
    elif len(result) == 2:
        newname=os.path.split(result[0])[-1].replace(
            'agilentg4502a_07_1__unc_edu__Level_3__unc_lowess_normalization_gene_level__data',
            'agilentg4502a_07__unc_edu__Level_3__unc_lowess_normalization_gene_level__data')
        merge_files(result, newname)
        return newname
    else:
        raise ValueError('extract_and_merge_files can only handle two files')
    
def format_firehose_rsem(filename, output):
    HYB_REF = "Hybridization REF"
    GENE_ID = "gene_id"
    DATA = arrayio.read(filename)
    assert DATA._row_order == [HYB_REF]
    assert DATA._col_order == ["_SAMPLE_NAME", GENE_ID]
    genes = DATA.row_names(HYB_REF)
    gene_symbols = [None] * len(genes)
    gene_ids = [None] * len(genes)
    for i in range(len(genes)):
        x = genes[i].split("|")
        assert len(x) == 2
        gene_symbol, gene_id = x
        if gene_symbol == "?":
            gene_symbol = ""
        gene_ids[i] = gene_id
        gene_symbols[i] = gene_symbol
    f = file(output,'w')
    header = ["Gene ID", "Gene Symbol"] + DATA.col_names("_SAMPLE_NAME")
    f.write("\t".join(header)+'\n')
    for i in range(DATA.nrow()):
        x = [gene_ids[i], gene_symbols[i]] + DATA._X[i]
        assert len(x) == len(header)
        f.write("\t".join(map(str, x))+'\n')
    f.close()


def format_firehose_exonexp(filename, output):
    HYB_REF = "Hybridization REF"
    GENE_ID = "exon"
    DATA = arrayio.read(filename)
    assert DATA._row_order == [HYB_REF]
    assert DATA._col_order == ["_SAMPLE_NAME", GENE_ID]
    # Line 2 of file is:
    # exon raw_counts median_length_normalized RPKM [...]
    # Only want RPKM columns.
    col_headers = DATA.col_names(GENE_ID)
    I_sig = []
    for i in range(len(col_headers)):
        if i % 3 == 2:
            assert col_headers[i] == "RPKM"
            I_sig.append(i)
        else:
            assert col_headers[i] != "RPKM"
    genes = DATA.row_names(HYB_REF)
    chroms = [None] * len(genes)
    starts = [None] * len(genes)
    ends = [None] * len(genes)
    strands = [None] * len(genes)
    for i in range(len(genes)):
        x = genes[i].split(":")
        assert len(x) == 3
        chrom, x, strand = x
        x = x.split("-")
        assert len(x) == 2
        start, end = x
        chroms[i] = chrom
        starts[i] = start
        ends[i] = end
        strands[i] = strand
    x = DATA.col_names("_SAMPLE_NAME")
    sample_names = [x[i] for i in I_sig]
    f = file(output,'w')
    header = ["Chrom", "Start", "End", "Strand"] + sample_names
    f.write("\t".join(header)+'\n')
    prev_line = None
    for i in range(DATA.nrow()):
        sig = [DATA._X[i][j] for j in I_sig]
        x = [chroms[i], starts[i], ends[i], strands[i]] + sig
        assert len(x) == len(header)
        if x == prev_line:
            continue
        prev_line = x
        f.write("\t".join(map(str, x))+'\n')
    f.close()

def annotate_firehose_methylation(filename, output):
    f=file(filename,'r')
    text = f.readlines(2)
    f.close()
    handle = text[1].split('\t')
    assert handle[:5] == ['Composite Element REF','Beta_value',
                           'Gene_Symbol','Chromosome','Genomic_Coordinate']
    f=file(filename,'r')
    all_symbols = {}
    symbols=[]
    for i in f.readlines():
        words = i.split('\t')
        symbol = words[2]
        symbols = symbol.split(";")
        for x in symbols:
            all_symbols[x] = 1
    f.close()
    all_symbols = sorted(all_symbols)
    #Look up all the symbols in the genefinder.
    symbol2id = {}
    genes = genefinder.find_many_genes(all_symbols, tax_id=9606)
    for (symbol, gene) in zip(all_symbols, genes):
        gene_id = gene[1]
        if gene_id is None:
            gene_id = ""
        symbol2id[symbol] = gene_id
    handle = filelib.read_row(filename, header=1)
    samples_names = [handle._header[i] for i in range(len(handle._header)) if not (i-1)%4]
    header = ["Probe.ID", "Gene.ID", "Gene.Symbol", "Chromosome",
              "Genomic.Coordinate"] + samples_names
    f = file(output, 'w')
    f.write("\t".join(header)+'\n')
    with open(filename) as FileObj:
        for lines in FileObj:
            if lines.startswith('Hybridization REF') or lines.startswith('Composite Element REF'):
                continue
            items = lines.split('\t')
            probe_id = items[0]
            Gene_symbols = items[2]
            Chormosome = items[3]
            Genomic_coordinate = items[4]
            values = [items[i] for i in range(len(items)) if not (i-1)%4]
            symbols = Gene_symbols.split(";")
            ids = [symbol2id.get(x, "") for x in symbols]
            #gene_symbol = ";".join(symbols)
            gene_id = ";".join(map(str, ids))
            row = [probe_id,gene_id,Gene_symbols,Chormosome,Genomic_coordinate]+values
            assert len(row) == len(header)
            f.write("\t".join(map(str, row))+'\n')
    f.close()


def format_firehose_mirna(filename, output):
    matrix = [x for x in filelib.read_cols(filename)]
    HYB_REF = "Hybridization REF"
    GENE_ID = "miRNA_ID"
    assert matrix
    assert matrix[0][0] == HYB_REF
    assert matrix[1][0] == GENE_ID
    header0 = matrix[0]
    header1 = matrix[1]
    for i in range(1, len(header1), 3):
        assert header1[i] == "read_count"
        assert header1[i+1] == "reads_per_million_miRNA_mapped"
        assert header1[i+2] == "cross-mapped"

    sample_name = [header0[i] for i in range(2, len(header0), 3)]
    header = ["miRNA ID"] + sample_name
    f = file(output, 'w')
    f.write("\t".join(header)+'\n')
    for i in range(2, len(matrix)):
        x = [matrix[i][j] for j in range(2, len(matrix[i]), 3)]
        x = [matrix[i][0]] + x
        assert len(x) == len(header)
        f.write("\t".join(x)+'\n')
    f.close()

def format_firehose_rppa(filename, output):
    COMP_REF = "Composite.Element.REF"
    COMP_REF_H = "Composite_Element_REF"
    iter = filelib.read_row(filename, header=1)
    assert iter._header[0] == COMP_REF
    f = file(output, 'w')
    header = ["Gene Symbol", "Gene ID", "Antibody"] + iter._header[1:]
    f.write("\t".join(header)+'\n')
    for d in iter:
        x = getattr(d, COMP_REF_H)
        x = x.split("|")
        assert len(x) == 2
        x, antibody = x
        gene_symbols = [x.strip() for x in x.split()]
        x = genefinder.find_many_genes(gene_symbols, tax_id="9606")
        x = [x[1] for x in x]
        x = [x for x in x if x]
        gene_ids = x
        if gene_symbols == ["CDC2"]:
            gene_ids = ["983"]
        
        gene_symbol_str = ";".join(gene_symbols)
        gene_id_str = ";".join(map(str, gene_ids))
        x = [gene_symbol_str, gene_id_str, antibody] + d._cols[1:]
        assert len(x) == len(header)
        f.write("\t".join(map(str, x))+'\n')
    f.close()
    
def format_firehose_gistic(filename, output):
    iter = filelib.read_row(filename, header=1)
    header = ["Gene ID", "Gene Symbol"] + iter._header[2:]
    f = file(output, 'w')
    f.write("\t".join(header))
    for d in iter:
        gene_symbol = d.Gene_Symbol
        gene_id = d.Locus_ID
        x = [gene_id, gene_symbol] + d._cols[2:]
        assert len(x) == len(header)
        f.write("\t".join(map(str, x))+'\n')
    f.close()

    
def format_rsem_isoforms(txt_file, outfile):
    from genomicode import arrayplatformlib
    
    M = arrayio.read(txt_file)
    # detect platform
    x = arrayplatformlib.score_matrix(M, min_score=0.8)
    assert x, "Cannot identify platform."
    header, platform = x.header, x.platform_name
    
    probe_ids = M.row_names(header)
    #if kg5, convert to kg7
    if platform == 'UCSC_human_hg19_kg5':
        new_platform = 'UCSC_human_hg19_kg7'
        kg7_ids = arrayannot.convert_probe_ids(probe_ids,new_platform)
        kg7_header = 'Hybridization REF kg7'
        M = make_matrix_new_ids(M, kg7_ids,kg7_header,1) 
        probe_ids = M.row_names(kg7_header)   
    # add LocusLink ids
    LocusLink_ids = arrayannot.convert_probe_ids(probe_ids,'Entrez_ID_human')
    gene_symbol_ids = arrayannot.convert_probe_ids(probe_ids,'Entrez_Symbol_human')
    newMatrix = make_matrix_new_ids(M,LocusLink_ids,'Entrez_ID_human',2)
    newMatrix = make_matrix_new_ids(newMatrix,gene_symbol_ids,'Entrez_Symbol_human',3)
    #get rid of scaled_estimate
    assert 'scaled_estimate' in newMatrix._col_names['isoform_id']
    assert 'raw_count' in newMatrix._col_names['isoform_id']
    col_names = {}
    col_names['_SAMPLE_NAME'] =[newMatrix._col_names['_SAMPLE_NAME'][i]
                                for i in range(len(newMatrix._col_names['_SAMPLE_NAME']))
                                if not i%2]
    row_names = newMatrix._row_names.copy()
    row_order = newMatrix._row_order[:]
    col_order = newMatrix._col_order[:]
    col_order.remove('isoform_id')
    synonyms = newMatrix._synonyms.copy()
    X = []
    for line in newMatrix._X:
        line = [line[i] for i in range(len(line)) if not i%2]
        X.append(line)
    x = Matrix.InMemoryMatrix(
            X, row_names=row_names, col_names=col_names,
            row_order=row_order, col_order=col_order, synonyms=synonyms)
    f = file(outfile,'w')
    arrayio.tab_delimited_format.write(x, f)
    f.close()

def format_affymetrix(filename, output):
    HYB_REF = "Hybridization REF"
    DATA = arrayio.read(filename)
    assert DATA._row_order == [HYB_REF]
    assert DATA._col_order == ["_SAMPLE_NAME", 'Composite Element REF']
    genes = DATA.row_names(HYB_REF)
    f = file(output,'w')
    header = [ "Gene Symbol"] + DATA.col_names("_SAMPLE_NAME")
    f.write("\t".join(header)+'\n')
    for i in range(DATA.nrow()):
        row = [j if j is not None else '' for j in DATA._X[i] ] 
        x = [genes[i]] + row
        assert len(x) == len(header)
        f.write("\t".join(map(str, x))+'\n')
    f.close()

def format_agilent(filename, output):
    HYB_REF = "Hybridization REF"
    DATA = arrayio.read(filename)
    assert DATA._row_order == [HYB_REF]
    assert DATA._col_order == ["_SAMPLE_NAME", 'Composite Element REF']
    genes = DATA.row_names(HYB_REF)
    f = file(output,'w')
    header = [ "Gene Symbol"] + DATA.col_names("_SAMPLE_NAME")
    f.write("\t".join(header)+'\n')
    for i in range(DATA.nrow()):
        row = [j if j is not None else '' for j in DATA._X[i] ] 
        x = [genes[i]] + row
        assert len(x) == len(header)
        f.write("\t".join(map(str, x))+'\n')
    f.close()

    
def make_matrix_new_ids(DATA,output_ids,header,index):
 # Make a matrix with the new IDs.
    X = DATA._X
    row_names = DATA._row_names.copy()
    row_order = DATA._row_order[:]
    col_names = DATA._col_names.copy()
    col_order = DATA._col_order[:]
    synonyms = DATA._synonyms.copy()
    row_order.insert(index,header)
    row_names[header] = output_ids
    # Write the outfile.
    x = Matrix.InMemoryMatrix(
        X, row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order, synonyms=synonyms)
    return x

def process_data(data, txt_file, outfile):
    if data == 'RSEM_genes':
        format_firehose_rsem(txt_file, outfile)
    elif data == 'RSEM_exons':
        format_firehose_exonexp(txt_file, outfile)
    elif data == 'humanmethylation450':
        annotate_firehose_methylation(txt_file, outfile)
    elif data == 'mirnaseq':
        format_firehose_mirna(txt_file, outfile)
    elif data == 'clinical':
        raise NotImplementedError("have not figured out how to process")
    elif data == 'rppa':
        format_firehose_rppa(txt_file, outfile)
    elif data == 'cnv_gistic2':
        format_firehose_gistic(txt_file, outfile)
    elif data == 'RSEM_isoforms':
        format_rsem_isoforms(txt_file, outfile)
    elif data == 'agilent':
        format_agilent(txt_file, outfile)
    elif data == 'affymetrix':
        format_affymetrix(txt_file,outfile)
    else:
        raise ValueError("the data type is not matched to our list")
    print 'processing finished '

    
def main():
    parser = argparse.ArgumentParser()
    
    parser.add_argument(
        '--disease', 
        help='Which disease to download data from.  Format: BRCA')
    parser.add_argument(
        '--date', 
        help='Download data from this date.  If not given, get the most '
        'recent.  Format: 20140715')
    x = sorted(datatype_match)
    x = ", ".join(x)
    parser.add_argument(
        '--data', 
        help="Which type of data to download.  Possibilities: %s" % x)
    
    parser.add_argument(
        'output', nargs="?", 
        help='output file for the processed data')
    
    parser.add_argument(
        '--list_dates', action='store_true',
        help='show the dates that have data.  If disease is given, show the '
        'dates for this disease')
    parser.add_argument(
        '--list_diseases', action='store_true',
        help='show the available diseases.  If date is not given, '
        'show for the latest date.')
    parser.add_argument(
        '--list_data', dest='list_data', action='store_const',
        const=True, default=False,
        help='show the data available for the specify disease and date, '
        'if date is not given, give the lastest one')

    parser.add_argument(
        '--download_only', action='store_true',
        help='Download the raw data file without processing.')
    parser.add_argument(
        '--download_and_extract', action='store_true',
        help='Download the raw data file with extracting.')
    parser.add_argument(
        '--process_only', action='store_true',
        help='Process a previously downloaded file.  Can be the original '
        '.tar.gz archive or unzip folder.')
    parser.add_argument('--input', help='input file for process')

    args = parser.parse_args()
    if args.data:
        assert args.data in datatype_match.keys(), \
               'the data you enter is not recognized'
    if args.date:
        assert re.search(r'[0-9]{8}', args.date)
    if args.download_only:
        assert args.disease, 'please specify the disease'
        assert args.data, 'please specify the data'
    if args.process_only:
        assert args.data, 'please specify the data'
    all_dates, all_dates_list = retrieve_all_dates()
    if args.list_dates:
        assert not args.date
        assert not args.data
        print "Dates"
        if args.disease:
            raise NotImplementedError
        else:
            for date in all_dates_list:
                print date
        return
    elif args.list_diseases:
        assert not args.disease
        date = all_dates_list[-1]
        if args.date:
            date = args.date
        all_diseases = retrieve_diseases(date)
        print "Diseases available at %s" % date
        if not args.data:
            for name in all_diseases:
                print name
            return
        else:
            for name in all_diseases:
                all_data = get_data_type(name, date)
                if args.data in all_data:
                    print name
            return
    elif args.list_data:
        assert args.disease, "disease must be given."
        date = all_dates_list[-1]
        if args.date:
            date = args.date
        all_data = get_data_type(args.disease, date)
        for d in all_data:
            print d
        return


    if args.process_only:
        assert args.input
        assert os.path.exists(args.input), '%s does not exists' % args.input
        #if instance(args.input
        txt_file = extract_and_merge_files([args.input],datatype2resource[args.data])
        process_data(args.data, txt_file, args.output)
    elif args.download_only:
        assert args.disease, "disease must be given."
        assert args.data, "data must be given."
        date = sorted(all_dates_list)[-1]
        if args.date:
            date = args.date
        download_file(args.disease,date,args.data)
    elif args.download_and_extract:
        assert args.disease, "disease must be given."
        assert args.data, "data must be given."
        date = sorted(all_dates_list)[-1]
        if args.date:
            date = args.date
        filenames = download_file(args.disease,date,args.data)
        txt_file = extract_and_merge_files(filenames,datatype2resource[args.data])
    else:
        assert args.disease, "Please specify a disease to download."
        assert args.data, "data must be given."
        assert args.output, "Please specify output path."
        date = all_dates_list[-1]
        if args.date:
            date = args.date
        diseases_in_date = retrieve_diseases(date)
        assert args.disease in diseases_in_date, (
            'the diesease %s is not found in the date %s' %(
            args.disease,date))
            
##        assert args.data in all_data,('%s is not found in %s for %s' %(
##            args.data,require_date,args.disease))
        
        filenames = download_file(args.disease, date, args.data)
        txt_file = extract_and_merge_files(
            filenames, datatype2resource[args.data])
        process_data(args.data, txt_file, args.output)
            
if __name__ == '__main__':
    main()
