#gpr_module.py
import gzip
import os


def check_gpr(fileloc):
    """check if a file is a gpr format"""
    entry = ['Type', 'DateTime', 'Settings', 'GalFile', 'PixelSize',
             'Wavelengths', 'ImageFiles', 'NormalizationMethod',
             'NormalizationFactors', 'JpegImage', 'StdDev',
             'RatioFormulations', 'Barcode', 'BackgroundSubtraction',
             'ImageOrigin', 'JpegOrigin', 'Creator', 'Scanner',
             'FocusPosition', 'Temperature', 'LinesAveraged', 'Comment',
             'PMTGain', 'ScanPower', 'LaserPower', 'LaserOnTime', 'Filters',
             'ScanRegion', 'Supplier']
    if fileloc.endswith('gz'):
        fileObj = gzip.GzipFile(fileloc, 'r')
    else:
        fileObj = file(fileloc, 'r')
    text = fileObj.readlines()
    fileObj.close()
    text = [i.replace('\"', '') for i in text]
    text = [i for i in text if not i == '']  # remove the empty line
    if not text[0].startswith('ATF'):  # 'the gpr requires a ATF number'
        return False
    a = text[1].split()
    # 'should be interger'
    if not (len(a) == 2 and a[0].isdigit() and a[1].isdigit()):
        return False
    i = -1
    startline = -1
    while i < len(text):
        i = i + 1
        a = text[i].split('\t')
        if len(a) > 2:
            startline = i
            break
    for line in text[2:startline]:
        line = line.split('=')
        if line[0] not in entry:  # '%s is not in gpr entry'%line[0]
            return False
    column_name = text[startline].split('\t')
    column_need = ['Name', 'ID', 'Block', 'Column', 'Row',
                   'Log Ratio (635/532)']
    for i in column_need:
        if i not in column_name:
            return False
    return True


def extract_gpr(fileloc, keep):
    filename = os.path.split(fileloc)[-1]
    if not filename.endswith('gpr.gz') and not filename.endswith('gpr'):
        raise ValueError('input file is not a gpr file')
    samplename = ''
    if filename.endswith('gpr.gz'):
        fileObj = gzip.GzipFile(fileloc, 'r')
        samplename = filename[:-7]
    elif filename.endswith('gpr'):
        fileObj = file(fileloc, 'r')
        samplename = filename.s[:-4]
    text = fileObj.readlines()
    fileObj.close()
    #get the index of content line of the gpr file
    i = -1
    startline = -1
    while i < len(text):
        i = i + 1
        a = text[i].split('\t')
        if len(a) > 2:
            startline = i
            break
    #find the column needed
    logratio = []
    column_name = text[startline].replace('\"', '')
    column_name = column_name.split('\t')
    a = column_name.index('Name')
    b = column_name.index('ID')
    c = column_name.index('Block')
    d = column_name.index('Column')
    k = column_name.index('Row')
    f = column_name.index('Log Ratio (635/532)')
    for i in range(startline, len(text)):
        h = text[i].replace('\"', '')
        g = h.split('\t')
        if g[f] == 'Error':
            x = ''
        else:
            x = g[f]
        if len(keep) <= i - startline:
            keep.append([g[b], g[a], g[c], g[d], g[k]])
        else:
            assert [g[b], g[a], g[c], g[d], g[k]] == keep[i - startline], (
                'row not match'
            )
        logratio.append(x)
    # replace the samplename to the "Log Ratio (635/532)"
    logratio[0] = samplename
    return logratio, keep


def extract_multiple(fileloc, keep):
    filename = os.path.split(fileloc)[-1]
    if not filename.endswith('gpr.gz') and not filename.endswith('gpr'):
        raise ValueError('input file is not a gpr file')
    if filename.endswith('gpr.gz'):
        fileObj = gzip.GzipFile(fileloc, 'r')
    elif filename.endswith('gpr'):
        fileObj = file(fileloc, 'r')
    text = fileObj.readlines()
    fileObj.close()
    #get the index of content line of the gpr file
    i = -1
    startline = -1
    while i < len(text):
        i = i + 1
        a = text[i].split('\t')
        if len(a) > 2:
            startline = i
            break
    #find the column needed
    red_signal = []
    green_signal = []
    red_back = []
    green_back = []
    column_name = text[startline].replace('\"', '')
    column_name = column_name.split('\t')
    a = column_name.index('Name')
    b = column_name.index('ID')
    c = column_name.index('Block')
    d = column_name.index('Column')
    k = column_name.index('Row')
    g1 = column_name.index('F635 Mean')
    g2 = column_name.index('F532 Mean')
    g3 = column_name.index('B635 Median')
    g4 = column_name.index('B532 Median')
    for i in range(startline, len(text)):
        h = text[i].replace('\"', '')
        g = h.split('\t')
        if len(keep) <= i - startline:
            keep.append([g[b], g[a], g[c], g[d], g[k]])
        else:
            assert [g[b], g[a], g[c], g[d], g[k]] == keep[i - startline], (
                'row not match'
            )
        red_signal.append(g[g1])
        green_signal.append(g[g2])
        red_back.append(g[g3])
        green_back.append(g[g4])
    return red_signal[1:], green_signal[1:], red_back[1:], green_back[1:], keep
