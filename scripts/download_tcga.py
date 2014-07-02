import urllib2
import re
import argparse
import os
from genomicode import parselib, filelib
import arrayio
from genomicode import genefinder,timer

datatype_match = {'RSEM_genes':'RSEM_genes_normalized__data.Level_3',
                  'RSEM_exons':'exon_expression__data.Level_3',
                  'humanmethylation450':
                  'Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3',
                  'mirnaseq':
                  'Merge_mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.Level_3',
                  'clinical':'Merge_Clinical.Level_1',
                  'rppa':'.RPPA_AnnotateWithGene.Level_3'}

def read_url(url):
    response = urllib2.urlopen(url)
    html = response.read()
    timer.wait(2,'tcga')
    return html

def get_all_dates():
    url = 'http://gdac.broadinstitute.org/runs/info/stddata__runs_list.html'
    html = read_url(url)
    all_dates = {}
    link_tags = parselib.get_tags_and_contents(html,'a')
    for i in link_tags:
        url = parselib.get_href(i)
        date = re.search(r'[0-9]{4}_[0-9]{2}_[0-9]{2}', url)
        if date:
           date = date.group(0)
           all_dates[date] = url   
    return all_dates


def download_file(disease,date,datatype):
    link = 'http://gdac.broadinstitute.org/runs/stddata__'+date+'/data/'+disease+'/'+date.replace('_','')+'/'
    newlinks = get_all_datas_on_page(link)
    for newlink in newlinks:
        if datatype_match[datatype] in newlink and disease+'-FFPE' not in newlink:
            data = read_url(link+newlink)
            with open(newlink, "wb") as code:
                 code.write(data)
            print 'finished download %s' %newlink
            return newlink
    assert ValueError('download fails')

    
def get_data_type(disease,date):
    link = 'http://gdac.broadinstitute.org/runs/stddata__'+date+'/data/'+disease+'/'+date.replace('_','')+'/'
    newlinks = get_all_datas_on_page(link)
    result = []
    for datatype in datatype_match:
        for newlink in newlinks:
            if datatype_match[datatype] in newlink and disease+'-FFPE' not in newlink:
                  if datatype not in result:
                      result.append(datatype)
    return result
                
            
def get_disease_lastest(urls):
    disease = []
    pattern = re.compile(r'http://gdac.broadinstitute.org/runs/stddata__[0-9_]*/[A-Z]*.html')
    lastest = None
    for url in urls:
        disease_link = re.findall(pattern, url)
        if not disease_link:
            continue
        for link in disease_link:
            result = re.search(r'[A-Z]+', link)
            lastest = re.search(r'[0-9]{4}_[0-9]{2}_[0-9]{2}',link).group(0)
            if result and result.group(0) not in disease:
               disease.append(result.group(0))
    return disease,lastest 


def get_all_datas_on_page(page):
    html = read_url(page)
    links = []
    link_tags = parselib.get_tags_and_contents(html,'a')
    for i in link_tags:
        url = parselib.get_href(i)
        if url.startswith("gdac.broadinstitute.org_") and url.endswith(".gz"):
                links.append(url)
    return links



def read_date_page(page):
    html = read_url(page)
    links = []
    link_tags = parselib.get_tags_and_contents(html,'a')
    for i in link_tags:
        url = parselib.get_href(i)
        if url and url.startswith('http'):
            links.append(url)
    return links


def extract_files(gzfile):
    assert gzfile.endswith('tar.gz')
    import tarfile
    tfile = tarfile.open(gzfile, 'r:gz')
    gzname = os.path.split(gzfile)[-1]
    newdir = os.path.join(os.getcwd(),gzname[:-7])
    tfile.extractall(newdir)
    folder = os.listdir(newdir)
    directory = os.path.join(newdir,folder[0])
    assert os.path.exists(directory)
    files = os.listdir(directory)
    for filename in files:
        if filename.endswith('txt') and filename != 'MANIFEST.txt':
            return os.path.join(directory,filename)
    return None    

    
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
            gene_symbol = ";".join(symbols)
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
            
def process_data(data,txt_file,outfile):
    if data == 'RSEM_genes':
            format_firehose_rsem(txt_file, outfile)
    elif data == 'RSEM_exons':
        format_firehose_exonexp(txt_file, outfile)
    elif data == 'humanmethylation450':
        annotate_firehose_methylation(txt_file, outfile)
    elif data == 'mirnaseq':
        format_firehose_mirna(txt_file, outfile)
    elif data == 'clinical':
        raise
    elif data == 'rppa':
        format_firehose_rppa(txt_file, outfile)
    else:
        raise
    print 'processing finished '
    
def main():
    parser = argparse.ArgumentParser(description='download_tcga')
    parser.add_argument(
        '--download_only', dest='download', action='store_const',
        const=True, default=False, help='Download the raw data file')
    parser.add_argument(
        '--process_only', dest='process', action='store_const',
        const=True, default=False, help='process the gz data file')
    parser.add_argument(
        '--input', dest='input',default=None,type=str,
        help='input file for process')
    parser.add_argument(
        '--list_diseases', dest='list_disease', action='store_const',
        const=True,default=False, help='show the available disease types')
        
    parser.add_argument(
        '--disease', dest='disease', default=None, type=str,
        help='which disease to download data from')
    
    parser.add_argument(
        '--list_dates', dest='list_dates', action='store_const',
        const=True, default=False,
        help='show the dates that have data,if disease is given, show the dates for this disease')

    parser.add_argument(
        '--date', dest='date', default=None, type=str,
        help='download data from this date, if not given, get the most recent.')

    parser.add_argument(
        '--list_data', dest='list_data', action='store_const',
        const=True, default=False,
        help='show the data available for the specify disease and date, if not, give the lastest one')

    parser.add_argument(
        '--data', dest='data', default=None, type=str,
        help='give the data type')
    
    parser.add_argument(
        'output', nargs="?", default=None, type=str,
        help='output file for the processed data')
   
    args = parser.parse_args()
    if args.data:
        assert args.data in datatype_match.keys(),'the data you enter is not recognized'
    if args.list_data:
        assert args.disease,'please specify the disease'
    if args.download:
        assert args.disease,'please specify the disease'
        assert args.data,'please specify the data'
    if args.process:
        assert args.data,'please specify the data'
        
    homepage = 'https://confluence.broadinstitute.org/display/GDAC/Dashboard-Stddata'
    new_links = read_date_page(homepage)
    lastest_diseases, lastest_date = get_disease_lastest(new_links)
    all_dates_dict = get_all_dates()
    all_dates_dict[lastest_date] = homepage
    all_dates = all_dates_dict.keys()
    all_dates.sort(reverse=True)
    all_diseases = lastest_diseases
    disease_in_date_dict = {}
    for date in all_dates:
        urls = read_date_page(all_dates_dict[date])
        diseases_in_date, x = get_disease_lastest(urls)
        disease_in_date_dict[date] = diseases_in_date
        all_diseases = list(set(all_diseases).union(set(diseases_in_date)))
    if args.list_disease:
        for name in all_diseases:
            print name
    if args.disease:
        assert args.disease in all_diseases
    if args.list_dates:
        if not args.disease:
            for date in all_dates:
                print date
        else:
            date_for_disease = []
            for date in all_dates_dict:
                if args.disease in disease_in_date_dict[date]:
                    date_for_disease.append(date)
            date_for_disease.sort(reverse=True)
            for date in date_for_disease:
                print date
    if args.date:
        assert re.search(r'[0-9]{4}_[0-9]{2}_[0-9]{2}',args.date)
        assert args.date in all_dates,'there is no data available for the date you enter'
        require_date = args.date
    else:
        require_date = lastest_date
    if args.process:
        assert args.input
        assert os.path.exists(args.input),'%s does not exists'%args.input
        txt_file = extract_files(args.input)
        process_data(args.data, txt_file, args.output)
        return
    all_data = get_data_type(args.disease,require_date)
    if args.list_data:
        for i in all_data:
            print i
    if args.download:
        download_file(args.disease,require_date,args.data)
    else:
        if args.list_disease or args.list_dates or args.list_data:
            return
        assert args.disease in disease_in_date_dict[require_date], (
            'the diesease %s is not found in the date %s' %(
                args.disease,args.require_date))
        assert args.data in all_data,('%s is not found in %s for %s' %(
            args.data,require_date,args.disease))
        download_filename = download_file(args.disease,require_date,args.data)
        txt_file = extract_files(download_filename)
        process_data(args.data, txt_file, args.output)
            
if __name__ == '__main__':
    main()
