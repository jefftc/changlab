#!/usr/bin/env python

def find_Rscript():
    import os
    
    search_paths = [
        "/opt/local/bin/Rscript",
        "/usr/local/bin/Rscript",
        "/usr/bin/Rscript",
        ]
    path = None
    for spath in search_paths:
        if spath is None or not os.path.exists(spath):
            continue
        path = spath
        break
    assert path is not None, "I could not find: Rscript"
    return path
    

def run_R(script):
    # Run an R script in the background.  Makes sure nothing prints to
    # the screen.
    import os
    import tempfile
    import subprocess
    from genomicode.dwdnorm import _safe_unlink

    # Bug: assumes that Rscript can be executed.
    assert type(script) is type("")

    rscript = find_Rscript()
    temp_path = "."
    script_file = None
    try:
        x, script_file = tempfile.mkstemp(dir=temp_path); os.close(x)
        open(script_file, 'w').write(script)

        command = [rscript, script_file]
        p = subprocess.Popen(
            command, bufsize=0, stdin=subprocess.PIPE,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
        w, r = p.stdin, p.stdout
        w.close()
        output = r.read()
    finally:
        _safe_unlink(script_file)
    return output



def parse_indexes(MATRIX, indexes, indexes_include_headers):
    from genomicode import parselib
    
    max_index = MATRIX.ncol()
    num_headers = len(MATRIX._row_names)
    assert max_index, "empty matrix"

    I = []
    for s, e in parselib.parse_ranges(indexes):
        if indexes_include_headers:
            s, e = s-num_headers, e-num_headers
        assert s >= 1, "Index out of range: %s" % s
        assert e <= max_index, "Index out of range: %s" % e
        s, e = s - 1, min(e, max_index)
        I.extend(range(s, e))
    return I


def make_sample_names(MATRIX):
    num = MATRIX.ncol()
    names = ['Sample%d' % (i+1) for i in range(num)]
    return names


def make_array_names(MATRIX):
    num = MATRIX.ncol()
    names = ['Array%d' % (i+1) for i in range(num)]
    return names


def write_EIF_file(MATRIX, filename):
    # write EIF file
    array_names = make_array_names(MATRIX)
    
    handle = open(filename, 'w')
    print >>handle, "\t".join(array_names)
    for row in MATRIX._X:
        for i in range(len(row)):
            if row[i] < 1:
                row[i] = 1
        print >>handle, "\t".join(map(str, row))


def write_SIF_file(MATRIX, filename, batches):
    sample_names = make_sample_names(MATRIX)
    array_names = make_array_names(MATRIX)

    assert len(sample_names) == len(array_names)
    assert len(sample_names) == len(batches)
    
    handle = open(filename, 'w')
    header = 'Array name', 'Sample name', 'Batch'
    print >>handle, "\t".join(header)
    for x in zip(array_names, sample_names, batches):
        print >>handle, "\t".join(map(str, x))


def main():
    import os
    import sys
    import argparse
    import tempfile
    import shutil
    
    import arrayio
    from Betsy import module_utils
    from genomicode import config

    usage = "combatnorm.py [options] expression_file"
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("expression_file", help="Matrix to normalize.")
    #parser.add_argument("normalized_file", help="Normalized matrix.")
    parser.add_argument(
        "--indexes1", 
        help="Which columns in batch 1, E.g. 1-5,8 (1-based, "
        "inclusive).  Columns leftover will be in the remaining batches.")
    parser.add_argument(
        "--indexes2",
        help="If three or more batches, then specify the samples for the "
        "second batch.")
    # Allow up to 4 groups for now.
    parser.add_argument("--indexes3")
    parser.add_argument(
        "--indexes_include_headers", default=False, action="store_true",
        help="If not given (default), then column 1 is the first column "
        "with data.  If given, then column 1 is the very first column in "
        "the file, including the headers.")
    
    args = parser.parse_args()
    assert os.path.exists(args.expression_file), \
        "File not found: %s" % args.expression_file
    assert args.indexes1
    if args.indexes2:
        assert args.indexes1
    if args.indexes3:
        assert args.indexes2


    MATRIX = arrayio.read(args.expression_file)

    # Figure out the batch for each sample.
    assert MATRIX.ncol(), "empty matrix"
    # Should be 1 to <num_groups>.
    batches = [None] * MATRIX.ncol()
    all_indexes = [args.indexes1, args.indexes2, args.indexes3]
    for i, indexes in enumerate(all_indexes):
        if not indexes:
            continue
        I = parse_indexes(MATRIX, indexes, args.indexes_include_headers)
        for j in I:
            assert j < len(batches)
            assert batches[j] is None, "overlapping: %d" % j
            batches[j] = i+1
    # The remaining samples are in the next batch.
    x = batches
    x = [x for x in batches if x]
    max_batch = max(x)
    for i in range(len(batches)):
        if batches[i] == None:
            batches[i] = max_batch + 1


    temp_path = None
    cwd = os.getcwd()
    try:
        temp_path = tempfile.mkdtemp(dir=".")
        os.chdir(temp_path)
        
        EIF_file = "EIF.dat"
        SIF_file = "SIF.dat"

        # Write the input files.
        write_EIF_file(MATRIX, EIF_file)
        write_SIF_file(MATRIX, SIF_file, batches)
    
        # Run combat.
        ComBat_R = config.combat
        assert os.path.exists(ComBat_R), 'File not found: %s' % ComBat_R
        
        assert os.path.exists(EIF_file)
        assert os.path.exists(SIF_file)
        x = [
            'source("%s")' % ComBat_R,
            'ComBat("%s", "%s", write=T)' % (EIF_file, SIF_file),
            ]
        script = "\n".join(x)
        run_R(script)
        #R = jmath.start_R()
        #R('source("%s")' % ComBat_R)
        #R('ComBat("%s", "%s", write=T)' % (EIF_file, SIF_file))
        outfile = "Adjusted_EIF.dat_.xls"
        assert module_utils.exists_nz(outfile), "File not found: %s" % outfile

        x = open(outfile).readlines()
        x = [x.rstrip("\r\n") for x in x]
        x = [x for x in x if x]
        x = [x.split("\t") for x in x]
        #header = x[0]
        norm = x[1:]
        MATRIX_norm = MATRIX.matrix()
        assert len(norm) == MATRIX_norm.nrow()
        for x in norm:
            assert len(x) == MATRIX_norm.ncol()
        MATRIX_norm._X = norm
    finally:
        if temp_path is not None and os.path.exists(temp_path):
            shutil.rmtree(temp_path)
        os.chdir(cwd)

    arrayio.tab_delimited_format.write(MATRIX_norm, sys.stdout)

 
if __name__=='__main__':
    main()
