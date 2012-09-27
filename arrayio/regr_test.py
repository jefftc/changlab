def print_fn(obj):
    return str(obj)

def test_tab_delimited_format():
    import StringIO
    from arrayio import tab_delimited_format

    filename = "samples/0159_cl.rma"
    gct_file = "samples/0988.gct"
    test(tab_delimited_format.is_format, (filename,), {}, "True")
    test(tab_delimited_format.is_format, (gct_file,), {}, "False")
    
    DATA = tab_delimited_format.read(open(filename))
    test(DATA.dim, (), {}, "(499, 19)")

    x = "['Probe.Set.ID', 'Description', 'LocusLink', 'Gene.Symbol']"
    test(print_fn, (DATA._row_order,), {}, x)

    # Test writing.  Need to read without conversion, so there will be
    # no differences due to floating point conversion.
    DATA = tab_delimited_format.read(open(filename), datatype=None)
    handle = StringIO.StringIO()
    tab_delimited_format.write(DATA, handle)
    test(handle.getvalue, (), {}, open(filename).read())

def test_pcl_format():
    import StringIO
    from arrayio import pcl_format

    filename = "samples/hallstrom.pcl"
    gct_file = "samples/0988.gct"
    test(pcl_format.is_format, (filename,), {}, "True")
    test(pcl_format.is_format, (gct_file,), {}, "False")
    
    DATA = pcl_format.read(open(filename))
    test(DATA.dim, (), {}, "(49, 67)")

    x = "['Gene.Symbol']"
    test(print_fn, (DATA._row_order,), {}, x)

    # Test writing.  Need to read without conversion, so there will be
    # no differences due to floating point conversion.
    DATA = pcl_format.read(open(filename), datatype=None)
    handle = StringIO.StringIO()
    pcl_format.write(DATA, handle)
    test(handle.getvalue, (), {}, open(filename).read())

def test_jeffs_format():
    import StringIO
    from arrayio import jeffs_format

    filename = "samples/0159_cl.rma"
    gct_file = "samples/0988.gct"
    test(jeffs_format.is_format, (filename,), {}, "True")
    test(jeffs_format.is_format, (gct_file,), {}, "False")

    DATA = jeffs_format.read(open(filename))
    test(DATA.dim, (), {}, "(499, 19)")

    x = "['Probe.Set.ID', 'Description', 'LocusLink', 'Gene.Symbol']"
    test(print_fn, (DATA._row_order,), {}, x)

    # Test writing.  Need to read without conversion, so there will be
    # no differences due to floating point conversion.
    DATA = jeffs_format.read(open(filename), datatype=None)
    handle = StringIO.StringIO()
    jeffs_format.write(DATA, handle)
    test(handle.getvalue, (), {}, open(filename).read())

def test_cdt_format():
    import StringIO
    from arrayio import cdt_format

    filename = "samples/hallstrom.cdt"
    gct_file = "samples/0988.gct"
    test(cdt_format.is_format, (filename,), {}, "True")
    test(cdt_format.is_format, (gct_file,), {}, "False")
    
    DATA = cdt_format.read(open(filename))
    test(DATA.dim, (), {}, "(49, 67)")

    x = "['GID', 'Gene.Symbol', 'NAME', 'GWEIGHT']"
    test(print_fn, (DATA._row_order,), {}, x)

    # Test writing.  Need to read without conversion, so there will be
    # no differences due to floating point conversion.
    DATA = cdt_format.read(open(filename), datatype=None)
    handle = StringIO.StringIO()
    cdt_format.write(DATA, handle)
    test(handle.getvalue, (), {}, open(filename).read())

def test(fn, args, keywds, standard):
    try:
        output = fn(*args, **keywds)
    except Exception, x:
        output = "exception"
        if output != standard:
            raise
        #raise
    status = "PASSED"
    if str(output) != str(standard):
        status = "FAILED"
        #print output
        #print standard
    MAXLEN = 20
    output = _hash(output)[:MAXLEN]
    standard = _hash(standard)[:MAXLEN]
    args = _hash(args)[:MAXLEN]
    x = status, args, output, standard
    print "\t".join(map(str, x))

def _hash(x):
    x = str(x)
    x = x.replace("\t", "\\t")
    x = x.replace("\n", "\\n")
    return x

def test_gct_format():
    import StringIO
    from arrayio import gct_format

    filename = "samples/0988.gct"
    file_ = "samples/0159_cl.rma"
    test(gct_format.is_format, (filename,), {}, "True")
    test(gct_format.is_format, (file_,), {}, "False")
    
    DATA = gct_format.read(open(filename))
    test(DATA.dim, (), {}, "(500, 59)")

    x = "['NAME', 'Description']"
    test(print_fn, (DATA._row_order,), {}, x)

    # Test writing.  Need to read without conversion, so there will be
    # no differences due to floating point conversion.
    DATA = gct_format.read(open(filename), datatype=None)
    handle = StringIO.StringIO()
    gct_format.write(DATA, handle)
    test(handle.getvalue, (), {}, open(filename).read())

def test_format_conversion():
    import StringIO
    import arrayio
    
    file_jeff = "samples/0159_cl.small.rma"
    file_pcl = "samples/0159_cl.small.pcl"
    file_gct = "samples/0159_cl.small.gct"

    # Test choose_format.
    fmt = arrayio.choose_format(file_jeff)
    test(print_fn, (fmt.__name__,), {}, "arrayio.jeffs_format")
    fmt = arrayio.choose_format(file_pcl)
    test(print_fn, (fmt.__name__,), {}, "arrayio.pcl_format")
    fmt = arrayio.choose_format(file_gct)
    test(print_fn, (fmt.__name__,), {}, "arrayio.gct_format")

    # Test guess_format.
    X_jeff = arrayio.read(file_jeff, datatype=None)
    X_pcl = arrayio.read(file_pcl, datatype=None)
    X_gct = arrayio.read(file_gct, datatype=None)
    fmt = arrayio.guess_format(X_jeff)
    test(print_fn, (fmt.__name__,), {}, "arrayio.jeffs_format")
    fmt = arrayio.guess_format(X_pcl)
    test(print_fn, (fmt.__name__,), {}, "arrayio.pcl_format")
    fmt = arrayio.guess_format(X_gct)
    test(print_fn, (fmt.__name__,), {}, "arrayio.gct_format")
    
    # Read the matrix.  No format conversion, or float conversion
    # might mess things up.
    X_jeff = arrayio.read(file_jeff, datatype=None)

    # Test convert.
    # _jeff_to_pcl
    handle = StringIO.StringIO()
    X_pcl = arrayio.convert(X_jeff, to_format=arrayio.pcl_format)
    arrayio.pcl_format.write(X_pcl, handle)
    #test(handle.getvalue, (), {}, open(file_pcl).read())
    
    #_jeff_to_gct
    handle = StringIO.StringIO()
    X_gct = arrayio.convert(X_jeff, to_format=arrayio.gct_format)
    arrayio.gct_format.write(X_gct, handle)
    #print handle.getvalue(),
    #test(handle.getvalue, (), {}, open(file_gct).read())

    #_gct_to_pcl
    handle = StringIO.StringIO()
    X_pcl = arrayio.convert(X_gct, to_format=arrayio.pcl_format)
    arrayio.pcl_format.write(X_pcl, handle)
    # _gct_to_pcl changes the gene id to "GeneID".  Fix this so that
    # we can compare against the gold standard file.  Everything else
    # should be the same.
    x = handle.getvalue()
    x = x.replace("GeneID", "Probe.Set.ID")
    test(print_fn, (x,), {}, open(file_pcl).read())
    
    #_pcl_to_gct
    handle = StringIO.StringIO()
    X_pcl = arrayio.read(file_pcl, datatype=None)
    X_gct = arrayio.convert(X_pcl, to_format=arrayio.gct_format)
    arrayio.gct_format.write(X_gct, handle)
    test(handle.getvalue, (), {}, open(file_gct).read())
    
    #_tdf_to_gct
    handle = StringIO.StringIO()
    X_pcl = arrayio.read(file_pcl, datatype=None)
    X_gct = arrayio.convert(
        X_pcl, from_format=arrayio.tab_delimited_format,
        to_format=arrayio.gct_format)
    arrayio.gct_format.write(X_gct, handle)
    test(handle.getvalue, (), {}, open(file_gct).read())
    

def main():
    test_tab_delimited_format()
    test_pcl_format()
    test_jeffs_format()
    test_cdt_format()
    test_gct_format()
    test_format_conversion()


if __name__ == '__main__':
    main()
