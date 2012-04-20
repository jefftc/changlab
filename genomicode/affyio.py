"""

Functions:
scan_cel
scan_celv3
scan_celv4
scan_celvcc1
scan_calvin_generic_data_file

guess_cel_version      Return a version number.
guess_cel_fn           Return a function.
extract_chip_name      Return the name of an array.

convert_cel_cc1_to_3

scan_bpmapv3

scan_cdf

"""

import os, sys

def scan_cel(filename):
    fn = guess_cel_fn(filename)
    return fn(filename)

def scan_celv3(filename):
    # Yields:
    # (SECTION, NAME, VALUE)
    # ("INTENSITY", "DATA", (X, Y, MEAN, STDEV, NPIXELS))
    # ("MASKS", "DATA", (X, Y))
    # ("OUTLIERS", "DATA", (X, Y))
    # ("MODIFIED", "DATA", (X, Y))
    import filelib
    
    assert type(filename) is type(""), "Need actual filename."
    handle = filelib.openfh(filename)  # in case of GZ file
    section = None
    for i, line in enumerate(handle):
        line = line.strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1]
        elif section == "INTENSITY" and line.find("=") < 0:
            x = line.strip().split()
            #if len(x) != 5:
            #    y = line.replace("\0", "")
            #    print repr(y)
            assert len(x) == 5, "Broken INTENSITY line: %s" % line.strip()
            x = int(x[0]), int(x[1]), float(x[2]), float(x[3]), int(x[4])
            yield section, "DATA", x
        elif section == "MASKS" and line.find("=") < 0:
            x = line.strip().split()
            assert len(x) == 2
            x = [int(x) for x in x]
            yield section, "DATA", x
        elif section == "OUTLIERS" and line.find("=") < 0:
            x = line.strip().split()
            assert len(x) == 2
            x = [int(x) for x in x]
            yield section, "DATA", x
        elif section == "MODIFIED" and line.find("=") < 0:
            x = line.strip().split()
            assert len(x) == 3
            x = int(x[0]), int(x[1]), float(x[2])
            yield section, "DATA", x
        else:
            assert section
            assert line.find("=") >= 0, line
            name, value = [x.strip() for x in line.split("=", 1)]
            yield section, name, value
            
    # If I opened this file, then close it.  gunzip might not die.
    handle.close()   

def scan_celv4(filename):
    # Yields:
    # (SECTION, NAME, VALUE)
    # ("INTENSITY", "DATA", (X, Y, MEAN, STDEV, NPIXELS))
    # ("MASKS", "DATA", (X, Y))
    # ("OUTLIERS", "DATA", (X, Y))
    # ("MODIFIED", "DATA", (X, Y))
    import struct
    import filelib
    
    # integer   32-bit signed integer
    # DWORD     32-bit unsigned integer
    # float     32-bit floating-point number
    # short     16-bit signed integer
    # little-endian
    def read(fmt):
        size = struct.calcsize(fmt)
        return struct.unpack(fmt, handle.read(size))

    assert type(filename) is type("")
    handle = filelib.openfh(filename, "rb")
    #handle.seek(0)

    magic, version = read("<ii")
    assert magic == 64
    assert version == 4
    yield "CEL", "Version", version

    num_cols, num_rows, num_cells = read("<iii")
    assert num_cells == num_cols * num_rows
    yield "HEADER", "Cols", num_cols
    yield "HEADER", "Rows", num_rows

    # The entire HEADER section of the CEL v3 files.
    length, = read("<i")
    header, = read("<%ds" % length)
    yield "HEADER", "Header", header

    length, = read("<i")
    algorithm, = read("<%ds" % length)
    yield "HEADER", "Algorithm", algorithm

    length, = read("<i")
    parameters, = read("<%ds" % length)
    yield "HEADER", "AlgorithmParameters", parameters

    cell_margin, num_outliers, num_masked, num_sub_grids = read("<iIIi")

    # Optimize the unpacking here.
    READ_SIZE = 100000   # 10 bytes each
    total_to_read = num_cells
    index = 0
    while total_to_read:
        n = min(total_to_read, READ_SIZE)
        fmt = "<" + "ffh"*n
        data = read(fmt)
        total_to_read -= n
        for i in range(0, len(data), 3):
            x = index % num_cols
            y = index / num_cols
            mean, stdev, npixels = data[i:i+3]
            yield "INTENSITY", "DATA", (x, y, mean, stdev, npixels)
            index += 1
    
    for i in range(num_masked):
        x, y = read("<hh")
        yield "MASKS", "DATA", (x, y)

    for i in range(num_outliers):
        x, y = read("<hh")
        yield "OUTLIERS", "DATA", (x, y)

    for i in range(num_sub_grids):
        row, col = read("<ii")
        x = read("<ffff")
        upper_left_x, upper_left_y, upper_right_x, upper_right_y = x
        x = read("<ffff")
        lower_left_x, lower_left_y, lower_right_x, lower_right_y = x
        x = read("<ffff")
        left_cell_pos, top_cell_pos, right_cell_pos, bottom_cell_pos = x

    #if type(filename) is type(""):
    handle.close()


def scan_celvcc1(filename):
    # Yields:
    # (SECTION, NAME, VALUE)
    # ("INTENSITY", "DATA", (X, Y, MEAN, STDEV, NPIXELS))
    # ("MASKS", "DATA", (X, Y))
    # ("OUTLIERS", "DATA", (X, Y))
    # ("MODIFIED", "DATA", (X, Y))
    #
    # Bug: Does not yield every piece of data from the file.  Will
    # need to manually add the ones that are of interest.

    # X is rows.
    # Y is columns.
    # X, Y coords start with (0, 0), (1, 0), (2, 0), etc.
    # Stored by row, then column.
    import filelib

    def scan_data_sets1(filename):
        # Yield list of tuples for each Data Set.
        data = []
        for x in scan_calvin_generic_data_file(filename):
            section, name, value = x
            if section != "DATA SET":
                continue
            if name == "NAME":
                if data:
                    yield data
                data = []
            data.append(x)
        if data:
            yield data
            
    def scan_data_sets2(filename):
        # Yield object with name, column_names, num_cols, num_rows, data.
        for data in scan_data_sets1(filename):
            dataset_name = None
            column_names = []
            values = []
            num_cols = num_rows = 0
            for x in data:
                section, name, value = x
                assert section == "DATA SET"
                if name == "TABLE":
                    values.append(value)
                elif name == "COLUMN":
                    column_names.append(value)
                elif name == "NAME":
                    dataset_name = value
                elif name == "NUM COLUMNS":
                    num_cols = value
                elif name == "NUM ROWS":
                    num_rows = value
            assert dataset_name is not None
            assert num_cols == len(column_names)
            assert num_cols*num_rows == len(values)

            x = filelib.GenericObject(
                name=dataset_name, column_names=column_names,
                num_cols=num_cols, num_rows=num_rows, data=values)
            yield x

    # Get the rows and cols on the chip.
    cel_rows = cel_cols = None
    for x in scan_calvin_generic_data_file(filename):
        section, name, value = x
        if section == "DATA GROUP":
            break
        if section != "DATA HEADER":
            continue
        yield section, name, value
        if name == "affymetrix-cel-rows":
            cel_rows = value
        elif name == "affymetrix-cel-cols":
            cel_cols = value
        #if cel_rows is not None and cel_cols is not None:
        #    break
    else:
        raise AssertionError, "Could not find chip dimensions"

    # Pull out the Data Sets.
    name2data = {}
    for x in scan_data_sets2(filename):
        name2data[x.name] = x

    assert "Intensity" in name2data
    assert "StdDev" in name2data
    assert "Pixel" in name2data
    x = len(name2data["Intensity"].data)
    assert x == name2data["Intensity"].num_rows
    assert x == len(name2data["StdDev"].data)
    assert x == len(name2data["Pixel"].data)
    assert x == cel_rows * cel_cols
    z = 0
    for y in range(cel_cols):
        for x in range(cel_rows):
            mean = name2data["Intensity"].data[z]
            stddev = name2data["StdDev"].data[z]
            npixels = name2data["Pixel"].data[z]
            z += 1
            yield "INTENSITY", "DATA", (x, y, mean, stddev, npixels)

    for i in range(name2data["Mask"].num_rows):
        data = name2data["Mask"].data
        x, y = data[i*2], data[i*2+1]
        yield "MASKS", "DATA", (x, y)

    for i in range(name2data["Outlier"].num_rows):
        data = name2data["Outlier"].data
        x, y = data[i*2], data[i*2+1]
        yield "OUTLIERS", "DATA", (x, y)

def scan_calvin_generic_data_file(filename):
    VALUE2TYPE = [
        "BYTE", "UBYTE", "SHORT", "USHORT", "INT", "UINT", "FLOAT",
        "STRING", "WSTRING"]
        
    def _read_raw(fmt):
        # fmt is a struct format string.
        import struct
        size = struct.calcsize(fmt)
        return struct.unpack(fmt, handle.read(size))
    
    def _read(type):
        # big-endian format
        if type == "BYTE":
            x, = _read_raw(">b")
        elif type == "UBYTE":
            x, = _read_raw(">B")
        elif type == "SHORT":
            x, = _read_raw(">h")
        elif type == "USHORT":
            x, = _read_raw(">H")
        elif type == "INT":
            x, = _read_raw(">i")
        elif type == "UINT":
            x, = _read_raw(">I")
        elif type == "FLOAT":
            x, = _read_raw(">f")
        elif type in ["STRING", "GUID", "VALUE"]:
            length, = _read_raw(">i")
            x, = _read_raw(">%ds" % length)
        elif type in ["WSTRING", "DATETIME", "LOCALE", "TYPE"]:
            length, = _read_raw(">i")
            x = handle.read(length*2)
            x = x.decode("utf-16BE").encode("utf-8")
        else:
            raise AssertionError, "Unknown type: %s" % type
        #print x
        return x

    def decode_mime(text, mime_type):
        import struct
        
        if mime_type == "text/ascii":
            # Strip NULL characters.
            x = [x for x in text if x != "\x00"]
            text = "".join(x)
        elif mime_type == "text/plain":
            x = text.decode("utf-16BE").encode("utf-8")
            x = [x for x in text if x != "\x00"]
            text = "".join(x)
        elif mime_type == "text/x-calvin-float":
            text, = struct.unpack(">f", text[:struct.calcsize("f")])
        elif mime_type == "text/x-calvin-integer-32":
            text, = struct.unpack(">i", text[:struct.calcsize("i")])
        elif mime_type == "text/x-calvin-unsigned-integer-16":
            text, = struct.unpack(">h", text[:struct.calcsize("h")])
        elif mime_type == "text/x-calvin-unsigned-integer-8":
            text, = struct.unpack(">h", "\x00" + text[0])
        else:
            raise AssertionError, "unhandled MIME type: %s" % mime_type
        return text

    # SECTION: File Header
    assert type(filename) is type("")
    # Do not accept .gz files because we need to seek to specific
    # locations.
    assert not filename.lower().endswith(".gz"), \
           "I cannot handle compressed files"
    handle = open(filename, "rb")
    #handle = filelib.openfh(filename, "rb")
    #handle.seek(0)
    magic = _read("UBYTE")
    version = _read("UBYTE")
    num_data_groups = _read("INT")
    pos_first_data_group = _read("UINT")
    assert magic == 59, repr(magic)
    assert version == 1

    yield "FILE HEADER", "NUM DATA GROUPS", num_data_groups

    # SECTION: Data Header
    # Read the header and all its parents.  Does not preserve the
    # lineage of the parents.
    num_data_headers = 1
    while num_data_headers:
        data_type_identifier = _read("GUID")
        unique_file_identifier = _read("GUID")
        yield "DATA HEADER", "DATA TYPE IDENTIFIER", data_type_identifier
        yield "DATA HEADER", "UNIQUE FILE IDENTIFIER", unique_file_identifier

        file_creation_date = _read("DATETIME")
        locale = _read("LOCALE")

        num_parameters = _read("INT")
        yield "DATA HEADER", "NUM PARAMETERS", num_parameters
        for i in range(num_parameters):
            name = _read("WSTRING")
            value = _read("VALUE")
            mime_type = _read("TYPE")
            value = decode_mime(value, mime_type)
            yield "DATA HEADER", name, value

        num_parent_file_headers = _read("INT")
        num_data_headers += num_parent_file_headers   # add the parents
        num_data_headers -= 1                         # I finished reading one.

    # SECTION: Data Group
    handle.seek(pos_first_data_group)
    for i in range(num_data_groups):
        pos_next_group = _read("UINT")
        pos_first_data_set = _read("UINT")
        num_data_sets = _read("INT")

        data_group_name = _read("WSTRING")
        yield "DATA GROUP", "NAME", data_group_name

        handle.seek(pos_first_data_set)
        for i in range(num_data_sets):
            # SECTION: Data Set
            pos_first_data_element = _read("UINT")
            pos_next_data_set = _read("UINT")

            data_set_name = _read("WSTRING")
            yield "DATA SET", "NAME", data_set_name

            num_parameters = _read("INT")
            yield "DATA SET", "NUM PARAMETERS", num_parameters
            for i in range(num_parameters):
                name = _read("WSTRING")
                value = _read("VALUE")
                mime_type = _read("TYPE")
                value = decode_mime(value, mime_type)
                yield "DATA SET", name, value

            num_columns = _read("UINT")
            yield "DATA SET", "NUM COLUMNS", num_columns
            columns = []
            for i in range(num_columns):
                col_name = _read("WSTRING")
                value_type = _read("BYTE")
                type_size = _read("INT")
                x = col_name, value_type, type_size
                columns.append(x)
                yield "DATA SET", "COLUMN", col_name

            num_rows = _read("UINT")
            yield "DATA SET", "NUM ROWS", num_rows
            for i in range(num_rows):
                for x in columns:
                    col_name, value_type, type_size = x
                    value = _read(VALUE2TYPE[value_type])
                    yield "DATA SET", "TABLE", value

            # END: Data Set
            handle.seek(pos_next_data_set)

        if pos_next_group != 0:
            handle.seek(pos_next_group)

    #if type(filename) is type(""):
    handle.close()

def scan_bpmapv3(filename):
    # Yields:
    # (SECTION, NAME, VALUE)
    # ("DESCRIPTION", "SEQUENCE_ID", <data>)
    # ("DESCRIPTION", "NAME", <data>)
    # ("DESCRIPTION", "TYPE", <data>)             0 (PM/MM); 1 (PM-only)
    # ("DESCRIPTION", "OFFSET", <data>)           file offset of POSITION_INFO
    # ("DESCRIPTION", "GROUP", <data>)
    # ("DESCRIPTION", "VERSION", <data>)
    # ("DESCRIPTION", "PARAMETER", (<name>, <value>))
    # ("POSITION_INFO", "SEQUENCE_ID", <data>)
    # ("POSITION_INFO", "PM_COORD", (<x>, <y>))   0-based coordinate
    # ("POSITION_INFO", "MM_COORD", (<x>, <y>))
    # ("POSITION_INFO", "PROBE_SEQ", <data>)
    # ("POSITION_INFO", "MATCH_SCORE", <data>)    always 1
    # ("POSITION_INFO", "PROBE_POS", <data>)      0-based position
    # ("POSITION_INFO", "STRAND", <data>)         1 target on +, 0 target on -
    import struct
    import filelib

    def read(fmt):
        size = struct.calcsize(fmt)
        return struct.unpack(fmt, handle.read(size))
    def read_string():
        length, = read(">I")
        return read(">%ds" % length)[0]

    # big-endian
    assert type(filename) is type("")
    handle = filelib.openfh(filename, "rb")
    #handle.seek(0)

    magic, = read(">8s")
    assert magic == "PHT7\r\n\x1a\n"

    # Bug in format: sometimes not stored as big-endian float.
    s = handle.read(4)
    version, = struct.unpack(">f", s)
    if int(version) not in [1, 2, 3]:
        version, = struct.unpack("<f", s)
    assert int(version) in [1, 2, 3]
    assert version == 3

    num_sequences, = read(">I")

    # SEQUENCE DESCRIPTION
    seq2nprobes = [None] * num_sequences
    seq2types = [None] * num_sequences
    for i in range(num_sequences):
        yield "DESCRIPTION", "SEQUENCE_ID", i
        sequence_name = read_string()
        yield "DESCRIPTION", "NAME", sequence_name

        probe_type, offset, num_probes = read(">III")
        seq2types[i] = probe_type
        seq2nprobes[i] = num_probes
        yield "DESCRIPTION", "TYPE", probe_type
        yield "DESCRIPTION", "OFFSET", offset
        group_name = read_string()
        yield "DESCRIPTION", "GROUP", group_name
        version = read_string()
        yield "DESCRIPTION", "VERSION", version

        num_params, = read(">I")
        for j in range(num_params):
            name = read_string()
            value = read_string()
            yield "DESCRIPTION", "PARAMETER", (name, value)

    # SEQUENCES
    for i in range(num_sequences):
        sequence_id, = read(">I")
        yield "POSITION_INFO", "SEQUENCE_ID", sequence_id
        for j in range(seq2nprobes[i]):
            pm_x, pm_y = read(">II")
            yield "POSITION_INFO", "PM_COORD", (pm_x, pm_y)
            if seq2types[i] == 0:
                mm_x, mm_y = read(">II")
                yield "POSITION_INFO", "MM_COORD", (mm_x, mm_y)

            probe_length, = read(">B")
            seq_code = read(">7B")
            x = read(">fIB")
            match_score, probe_pos, strand = x

            # Convert the code to a sequence.
            seq_str = []
            for k in range(len(seq_code)):
                s = "ACGT"
                seq_str.append(s[seq_code[k] >> 6 & 3])
                seq_str.append(s[seq_code[k] >> 4 & 3])
                seq_str.append(s[seq_code[k] >> 2 & 3])
                seq_str.append(s[seq_code[k] >> 0 & 3])
            seq_str = seq_str[:probe_length]
            seq_str = "".join(seq_str)
            yield "POSITION_INFO", "PROBE_SEQ", seq_str
            yield "POSITION_INFO", "MATCH_SCORE", match_score
            yield "POSITION_INFO", "PROBE_POS", probe_pos
            yield "POSITION_INFO", "STRAND", strand

    #if type(filename) is type(""):
    handle.close()
        
def guess_cel_version(filename):
    # Returns:
    # v3   Version 3 from MAS software.
    # v4   Version 4 from GCOS software.
    # cc1  Command Console version 1.
    import struct
    import filelib

    # Guess the version from the beginning of the file.

    # I need to be able to read from the start of the file.  If I
    # accept a file handle, it's not guaranteed to be at the start of
    # the file.  I can try to seek to the beginning of the file, but
    # this will fail for some files, e.g. gzip'd files.  It's easiest
    # just to not allow file handles.
    assert type(filename) is type("")
    handle = filelib.openfh(filename, "rb")
    #handle.seek(0)   # in case filename was a file handle
    data = handle.read(100)
    handle.close()   # close or gunzip may not die

    # Check to see if it has the magic numbers for version 4.
    size = struct.calcsize("<ii")
    magic, version = struct.unpack("<ii", data[:size])
    if magic == 64 and version == 4:
        return "v4"

    # Check to see if it has the magic numbers for Command Console
    # version 1.
    size = struct.calcsize(">BB")
    magic, version = struct.unpack(">BB", data[:size])
    if magic == 59 and version == 1:
        return "cc1"

    # See if it looks like version 3.
    # [CEL]
    # Version=3
    s = "[CEL]\nVersion=3"
    d = data
    d = d.replace("\r\n", "\n")
    d = d.replace("\r", "\n")
    if d[:len(s)] == s:
        return "v3"

    raise AssertionError, "Unable to guess CEL version for file %s" % filename

def convert_cel_cc1_to_3(filename, outhandle=None):
    outhandle = outhandle or sys.stdout

    #data = [x for x in scan_calvin_generic_data_file(filename)]
    data = []
    for x in scan_calvin_generic_data_file(filename):
        section, name, value = x
        #if section not in ["FILE HEADER", "DATA HEADER"]:
        #    break
        data.append(x)

    def getvalue(section, name, default=None):
        for s, n, v in data:
            if s == section and n == name:
                return v
        if default is not None:
            return default
        raise AssertionError, "not found: %s %s" % (section, name)

    def getdataset(name):
        in_dataset = 0
        ds = []
        for x in data:
            s, n, v = x
            if s != "DATA SET":
                continue
            if n == "NAME" and v != name:
                in_dataset = 0
                if len(ds):
                    break
            if n == "NAME" and v == name:
                assert not in_dataset
                in_dataset = 1
            if in_dataset:
                ds.append(x)
        assert ds, "not found: %s" % name
        return ds


    # CEL section.
    print >>outhandle, "[CEL]"
    print >>outhandle, "Version=3"
    print >>outhandle

    # HEADER section.
    cel_cols = getvalue("DATA HEADER", "affymetrix-cel-cols")
    cel_rows = getvalue("DATA HEADER", "affymetrix-cel-rows")
    print >>outhandle, "[HEADER]"
    print >>outhandle, "Cols=%d" % cel_cols
    print >>outhandle, "Rows=%d" % cel_rows
    print >>outhandle, "TotalX=%d" % cel_cols
    print >>outhandle, "TotalY=%d" % cel_rows
    print >>outhandle, "OffsetX=0"
    print >>outhandle, "OffsetY=0"
    x = getvalue("DATA HEADER", "affymetrix-algorithm-param-GridULX")
    y = getvalue("DATA HEADER", "affymetrix-algorithm-param-GridULY")
    print >>outhandle, "GridCornerUL=%d %d" % (x, y)
    x = getvalue("DATA HEADER", "affymetrix-algorithm-param-GridURX")
    y = getvalue("DATA HEADER", "affymetrix-algorithm-param-GridURY")
    print >>outhandle, "GridCornerUR=%d %d" % (x, y)
    x = getvalue("DATA HEADER", "affymetrix-algorithm-param-GridLRX")
    y = getvalue("DATA HEADER", "affymetrix-algorithm-param-GridLRY")
    print >>outhandle, "GridCornerLR=%d %d" % (x, y)
    x = getvalue("DATA HEADER", "affymetrix-algorithm-param-GridLLX")
    y = getvalue("DATA HEADER", "affymetrix-algorithm-param-GridLLY")
    print >>outhandle, "GridCornerLL=%d %d" % (x, y)
    print >>outhandle, "Axis-invertX=0"
    print >>outhandle, "AxisInvertY=0"
    print >>outhandle, "swapXY=0"

    # In CC1 files, the platform (e.g. HG-U133_Plus_2) can be stored
    # in several places.  Try them all.
    platform_names = [
        "affymetrix-dat-header", "affymetrix-partial-dat-header",
        "affymetrix-array-type", "affymetrix-created-arraytype",
        ]
    for name in platform_names:
        x = getvalue("DATA HEADER", name, "").strip()
        if not x:
            continue
        print >>outhandle, "DatHeader=%s" % x
        break
        
    print >>outhandle, "Algorithm=%s" % getvalue(
        "DATA HEADER", "affymetrix-algorithm-name")
    params = []
    for section, name, value in data:
        if section not in ["FILE HEADER", "DATA HEADER"]:
            break
        x = "affymetrix-algorithm-param-"
        if not name.startswith(x):
            continue
        n = name[len(x):]
        params.append("%s:%s" % (n, value))
    print >>outhandle, "AlgorithmParameters=%s" % ";".join(params)
    print >>outhandle
    
    ds_mean = getdataset("Intensity")
    ds_stdv = getdataset("StdDev")
    ds_npix = getdataset("Pixel")
    tab_mean = [x for x in ds_mean if x[1] == "TABLE"]
    tab_stdv = [x for x in ds_stdv if x[1] == "TABLE"]
    tab_npix = [x for x in ds_npix if x[1] == "TABLE"]
    assert len(tab_mean) == len(tab_stdv)
    assert len(tab_mean) == len(tab_npix)
    assert len(tab_mean) == cel_rows * cel_cols
    print >>outhandle, "[INTENSITY]"
    print >>outhandle, "NumberCells=%d" % len(tab_mean)
    print >>outhandle, "CellHeader=X\tY\tMEAN\tSTDV\tNPIXELS"
    z = 0
    for y in range(cel_cols):
        for x in range(cel_rows):
            x = x, y, tab_mean[z][2], tab_stdv[z][2], tab_npix[z][2]
            print >>outhandle, "%3d\t%3d\t%g\t%g\t%3d" % x
            z += 1
    print >>outhandle

    print >>outhandle, "[MASKS]"
    ds = getdataset("Mask")
    tab = [x for x in ds if x[1] == "TABLE"]
    print >>outhandle, "NumberCells=%d" % (len(tab)/2)
    print >>outhandle, "CellHeader=X\tY"
    for i in range(0, len(tab), 2):
        print >>outhandle, "%d\t%d" % (tab[i][2], tab[i+1][2])
    print >>outhandle

    print >>outhandle, "[OUTLIERS]"
    ds = getdataset("Outlier")
    tab = [x for x in ds if x[1] == "TABLE"]
    print >>outhandle, "NumberCells=%d" % (len(tab)/2)
    print >>outhandle, "CellHeader=X\tY"
    for i in range(0, len(tab), 2):
        print >>outhandle, "%d\t%d" % (tab[i][2], tab[i+1][2])
    print >>outhandle
    
    print >>outhandle, "[MODIFIED]"
    print >>outhandle, "NumberCells=0"
    print >>outhandle, "CellHeader=X\tY\tORIGMEAN"

def guess_cel_fn(filename):
    version = guess_cel_version(filename)
    if version == "v3":
        return scan_celv3
    elif version == "v4":
        return scan_celv4
    elif version == "cc1":
        return scan_celvcc1
    raise AssertionError, "no scanner for this version"

def extract_chip_name(filename):
    # Return the name of the chip or None.
    #
    # Example chip names:
    # HG-U133_Plus_2
    # HG-U133A
    # HG-U133A_2
    # HG-U133B
    # MOE430A
    # RAE230A
    # Mouse430A_2
    fn = guess_cel_fn(filename)

    chip_name = None
    for x in fn(filename):
        section, name, value = x

        # CEL v3, found in HEADER section:
        # DatHeader=[4..46107]  Yang 80 U133B:CLS=4733 RWS=4733 XIN=3  \
        # YIN=3  VE=17        2.0 07/10/03 15:42:44    ^T  ^T          \
        # HG-U133B.1sq ^T  ^T  ^T  ^T  ^T  ^T  ^T  ^T  ^T 6
        #
        # CEL vcc1, founder in DATA HEADER section:
        # affymetrix-array-type  HG-U133_Plus_2

        if section == "INTENSITY":
            break
        if section not in ["HEADER", "DATA HEADER"]:
            continue
        if name not in ["DatHeader", "Header", "affymetrix-array-type"]:
            continue

        # Raw chip names:
        # HG-U133_Plus_2.1sq
        # HG-U133A.1sq
        # HG-U133B.1sq
        # MOE430A.1sq
        # Mouse430A_2.1sq
        x = value.split()
        if len(x) > 1:
            x = [x for x in x if x.find("sq") >= 0]
        if not x:
            continue
        chip_name = x[0]
        # Remove the .1sq.
        chip_name = chip_name.replace(".1sq", "")
        break
    return chip_name

def scan_cdf(filename):
    # Yields:
    # (SECTION, NAME, VALUE)
    # ("CDF", "Version", VERSION)
    # ("Chip", "Name", NAME_OF_ARRAY)
    # ("Chip", "Rows", NUM_ROWS)
    # ("Chip", "Cols", NUM_COLS)
    #
    # ("Unit<x>", "UnitType", <type>)    <type> 3 is Expression.
    # ("Unit<x>", "NumAtoms", <num>)     Num probes, match/mismatch count as 1.
    # ("Unit<x>", "NumCells", <num>)     Num probes, match/mismatch count as 2.
    # ("Unit<x>", "NumberBlocks", <num>) Num blocks in probe set.
    # ("Unit_Block<x>", "CellHeader", <values>)  Header of table.
    #   X      X coordinate of cell.
    #   Y      Y coordinate of cell.
    #   QUAL   Probe set name.  For Genotyping units, includes allele.
    #   EXPOS  For Expression, ranges from [0, NumAtoms-1]
    #   ATOM   For Expression, same as EXPOS.  Groups together match/mismatch.
    #   INDEX  Used to look up cell data in CEL file.
    #   POS    Indexes within probe where mismatch occurs.
    #   PBASE  Base of probe at substitution position.
    #   TBASE  Base of target.  Is a PM probe if PBASE != TBASE.
    #          Otherwise, is a MM probe.  Seems flipped to me?  I guess
    #          it's a PM because it interrogates the complement.
    import filelib
    
    # ("Unit_Block<x>", "Cell<i>", <values>)  
    
    section = None
    Unit_Block_Cell_int_indexes = None
    assert type(filename) is type("")
    handle = filelib.openfh(filename)
    for i, line in enumerate(handle):
        line = line.strip("\r\n")
        if not line:
            continue
        
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1]
            continue

        assert section
        assert line.find("=") >= 0, line
        name, value = [x.strip() for x in line.split("=", 1)]

        # Do some cleaning up of the results.
        Chip_to_int = [
            "Rows", "Cols", "NumberOfUnits", "MaxUnit", "NumQCUnits"]
        Unit_to_int = [
            "Direction", "NumAtoms", "NumCells", "UnitNumber", "UnitType",
            "NumberBlocks"]
        Unit_Block_to_int = [
            "BlockNumber", "NumAtoms", "NumCells", "StartPosition",
            "StopPosition"]
        Unit_Block_Cell_to_int = [
            "X", "Y", "EXPOS", "POS", "ATOM", "INDEX", "CODONIND", "CODON",
            "REGIONTYPE"]
        
        if section == "Chip" and name in Chip_to_int:
            value = int(value)
        elif section.startswith("QC") and name == "CellHeader":
            value = tuple(value.split("\t"))
            assert len(value) == 6
        elif section.startswith("QC") and name.startswith("Cell"):
            value = value.split("\t")
            assert len(value) == 6
            x0, x1, x2, x3, x4, x5 = value
            value = int(x0), int(x1), x2, int(x3), int(x4), int(x5)
        elif(section.startswith("Unit") and section.find("Block") < 0 and
             name in Unit_to_int):
            value = int(value)
        elif(section.startswith("Unit") and section.find("Block") >= 0 and
             name in Unit_Block_to_int):
            value = int(value)
        elif(section.startswith("Unit") and section.find("Block") >= 0 and
             name == "CellHeader"):
            value = tuple(value.split("\t"))
            Unit_Block_Cell_int_indexes = [
                i for (i, x) in enumerate(value)
                if x in Unit_Block_Cell_to_int]
        elif(section.startswith("Unit") and section.find("Block") >= 0 and
             name.startswith("Cell")):
            value = value.split("\t")
            assert Unit_Block_Cell_int_indexes
            for i in Unit_Block_Cell_int_indexes:
                value[i] = int(value[i])
            value = tuple(value)

        yield section, name, value
    if type(filename) is type(""):
        handle.close()
