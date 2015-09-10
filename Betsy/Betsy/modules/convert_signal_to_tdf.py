from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        """check an input file is xls or xlsx format"""
        import arrayio
        
        in_filename = in_data.identifier
        # Why is this necessary?
        #try:
        #    x = userfile._unhash_storefile(in_data.identifier)
        #    real_name = x[1]
        #except:
        #    pass
        
        #if (in_data.identifier.endswith('.gz') or in_filename.endswith('.gz')):
        #    unzip_file = module_utils.gunzip(in_data.identifier)
        #else:
        #    unzip_file = in_data.identifier

        ## M = None
        ## xls_file = None
        ## txt_file = unzip_file
        ## try:
        ##     xlrd.open_workbook(unzip_file)
        ##     xls_file = 'tmp.xls'
        ## # XLRDError?  Is this a bug?  This is not the way to catch exception.
        ## except Exception, XLRDError:
        ##     try:
        ##         # Test this.  book not used?
        ##         book = openpyxl.load_workbook(unzip_file)
        ##         xls_file = 'tmp.xlsx'
        ##     except Exception, InvalidFileException:
        ##         xls_file = None
        ##     except (SystemError, MemoryError, KeyError), x:
        ##         raise
    
        
        ## if xls_file:
        ##     shutil.copyfile(unzip_file, xls_file)
        ##     xls2txt_path = config.xls2txt
        ##     xls2txt_BIN = module_utils.which(xls2txt_path)
        ##     assert xls2txt_BIN, 'cannot find the %s' % xls2txt_path
        ##     f = file('tmp1.txt', 'w')
        ##     command = ['python', xls2txt_BIN, xls_file]
        ##     process = subprocess.Popen(command,
        ##                                shell=False,
        ##                                stdout=f,
        ##                                stderr=subprocess.PIPE)
        ##     error_message = process.communicate()[1]
        ##     if error_message:
        ##         raise ValueError(error_message)
        ##     os.remove(xls_file)
        ##     f.close()
        ##     txt_file = 'tmp1.txt'

        to_format = arrayio.tdf
        MATRIX = arrayio.read(in_filename)
        MATRIX_c = arrayio.convert(MATRIX, to_format=to_format)
        to_format.write(MATRIX_c, open(outfile, 'w'))
        
        #M = guess_and_change_gct_header(txt_file)
        #M_c = arrayio.convert(M, to_format=arrayio.tab_delimited_format)
        #f = file(outfile, 'w')
        #arrayio.tab_delimited_format.write(M_c, f)
        #f.close()


    def name_outfile(self, antecedents, user_options):
        return "signal.tdf"
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #filename = 'signal_' + original_file + '.tdf'
        #return filename


## def guess_and_change_gct_header(filename):
##     import arrayio
##     from Betsy import module_utils
##     from genomicode import arrayplatformlib
    
##     CATEGORY2HEADER = {
##         arrayplatformlib.PROBE_ID : "Probe ID",
##         arrayplatformlib.GENE_ID : "Gene ID",
##         arrayplatformlib.GENE_SYMBOL : "Gene Symbol",
##         arrayplatformlib.DESCRIPTION : "Description",
##         }

##     M_name = arrayio.choose_format(filename)
##     M = arrayio.read(filename)
##     #ids = M._row_order
##     if not M_name.__name__ == 'arrayio.gct_format':
##         return M
    
##     all_platforms = arrayplatformlib.identify_all_platforms_of_matrix(M)
##     if not all_platforms:
##         return M

##     old_header, platform_name = all_platforms[0]

##     platform = arrayplatformlib.find_platform_by_name(platform_name)
##     assert platform
##     assert platform.category is not None, "No category: %s" % platform_name
##     assert platform.category in CATEGORY2HEADER, \
##            "Unknkown category: %s" % platform.category
##     new_header = CATEGORY2HEADER[platform.category]
##     #new_header = platform2header[platform]
##     M = module_utils.replace_matrix_header(M, old_header, new_header)
##     return M


## platform2header = {
##     'Agilent_Human1A': 'Probe ID',
##     'HG_U133A_2': 'Probe ID',
##     'HG_U133A': 'Probe ID',
##     'HG_U133B': 'Probe ID',
##     'HG_U133_Plus_2': 'Probe ID',
##     'HG_U95A': 'Probe ID',
##     'HG_U95Av2': 'Probe ID',
##     'Hu35KsubA': 'Probe ID',
##     'Hu35KsubB': 'Probe ID',
##     'Hu35KsubC': 'Probe ID',
##     'Hu35KsubD': 'Probe ID',
##     'Hu6800': 'Probe ID',
##     'HumanHT_12_control': 'Probe ID',
##     'HumanHT_12': 'Probe ID',
##     'MG_U74Av2': 'Probe ID',
##     'MG_U74Bv2': 'Probe ID',
##     'MG_U74Cv2': 'Probe ID',
##     'Mouse430_2': 'Probe ID',
##     'Mouse430A_2': 'Probe ID',
##     'MouseRef_8_control': 'Probe ID',
##     'MouseRef_8': 'Probe ID',
##     'Mu11KsubA': 'Probe ID',
##     'Mu11KsubB': 'Probe ID',
##     'RAE230A': 'Probe ID',
##     'RG_U34A': 'Probe ID',
##     'Entrez_ID_human': 'Gene ID',
##     'Entrez_ID_mouse': 'Gene ID',
##     'Entrez_symbol_human': 'Gene symbol',
##     'Entrez_symbol_mouse': 'Gene symbol'
## }
