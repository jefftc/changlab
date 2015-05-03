"""

Methods:
run
name_outfile
make_unique_hash
get_out_attributes
find_antecedents

"""


# data_node   Data object.
# parameters  For making a Data object.
# user_input  dictionary of parameters given by the user.

def run(data_node, parameters, user_input, network, num_cores):
    """given a database ID and GPLID, get the series matrix file"""
    import os
    import shutil

    from genomicode import geolib

    from Betsy import module_utils
    from Betsy import bie3
    from Betsy import rulebase
    
    GSEID = user_input['GSEID']
    GPLID = user_input.get("GPLID")
    assert GSEID.startswith('GSE'), 'GSEID %s is not correct' % GSEID
    assert not GPLID or GPLID.startswith('GPL'), \
           'GPLID %s is not correct' % GPLID
    
    outfile = name_outfile(data_node, user_input)
    outhandle = open(outfile, 'w')
    geolib.download_seriesmatrix_file(outhandle, GSEID, GPLID)
    outhandle.close()
    #if not os.path.exists(outfile):
    #    os.mkdir(outfile)
    #matrix_files = get_seriesmatrix_file(GSEID, GPLID)
    #for matrix_file in matrix_files:
    #    newmatrix_filename = os.path.split(matrix_file)[-1]
    #    shutil.copyfile(matrix_file, os.path.join(outfile, newmatrix_filename))
    #assert module_utils.exists_nz(outfile), (
    #    'the output file %s for download_geo_dseriesmatrix fails' % outfile
    #)
    
    out_node = bie3.Data(rulebase.GEOSeriesMatrixFile, **parameters)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


# Not sure exactly what this does.  Looks like it makes a unique name
# for the output of this module.
def name_outfile(data_node, user_input):
    import os
    from Betsy import module_utils

    # BUG: what if there are multiple GPLIDs?
    original_file = module_utils.get_inputid(user_input['GSEID'])
    filename = "download_geo_seriesmatrix_%s" % original_file
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def make_unique_hash(data_node, pipeline, parameters, user_input):
    GSEID = user_input["GSEID"]
    GPLID = user_input.get("GPLID", "")
    return "%s_%s" % (GSEID, GPLID)


def get_out_attributes(parameters, data_node):
    return parameters


def find_antecedents(network, module_id, data_nodes, parameters,
                     user_attributes):
    from Betsy import module_utils
    
    data_node = module_utils.get_identifier(network, module_id, data_nodes,
                                            user_attributes)
    return data_node


## def get_seriesmatrix_file(GSEID, GPLID=None):
##     'download series matrix and unzip'
##     from ftplib import FTP
##     import gzip
    
##     try:
##         ftp = FTP('ftp.ncbi.nih.gov')
##         ftp.login()
##     except Exception, e:
##         raise ValueError(e)
##     try:
##         ftp.cwd('pub/geo/DATA/SeriesMatrix/' + GSEID)
##     except FTP.error_perm, x:
##         if str(x).find('No such file') >= 0:
##             raise AssertionError('cannot find the %s' % path)
##     entry = []
##     ftp.retrlines('NLST', entry.append)
##     platform_txtfiles = []
##     for platform_filename in entry:
##         if GPLID and GPLID not in platform_filename:
##             continue
            
##         f = open(platform_filename, 'wb')
##         ftp.retrbinary('RETR ' + platform_filename, f.write)
##         f.close()
##         platform_txtfile = platform_filename[:-3]
##         assert not os.path.exists(platform_txtfile), (
##             'the seriesmatrix file %s already exists' % platform_txtfile
##         )
##         #unzip the gz data
##         fileObj = gzip.GzipFile(platform_filename, 'rb')
##         fileObjOut = file(platform_txtfile, 'wb')
##         while 1:
##             line = fileObj.readline()
##             if line == '':
##                 break
##             fileObjOut.write(line)
##         fileObj.close()
##         fileObjOut.close()
##         os.remove(platform_filename)
##         assert os.path.exists(platform_txtfile), (
##             'the unzip %s in download_geo_dataset_GPL fails' % platform_txtfile
##         )
##         platform_txtfiles.append(os.path.realpath(platform_txtfile))
##     ftp.close()
##     return platform_txtfiles
