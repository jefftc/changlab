"""

Functions:
score_tfbs
patser
parse

read_matrix
label_matrix
is_vertical_matrix
is_weight_matrix
is_labelled_matrix

"""
import os, sys

def score_tfbs(
    sequence, matrix_file, patser_bin=None, alphabet_file=None,
    num_jobs=None):
    # Return list of sequence_num, matrix_num, 0-based position,
    # strand, score, ln_p_value.  sequence can be a list of sequences.
    # matrix_file can be a single file or a list of matrix files.
    import tempfile
    import subprocess
    import config
    import parallel

    # In practice, parallelization of these jobs don't help things run
    # faster unless there are hundreds of matrices.
    patser_bin = patser_bin or config.patser
    alphabet_file = alphabet_file or config.patser_alphabet
    num_jobs = num_jobs or 1

    assert num_jobs >= 1

    sequences = sequence
    if type(sequence) is type(""):
        sequences = [sequence]
    matrix_files = matrix_file
    if type(matrix_file) is type(""):
        matrix_files = [matrix_file]
    for filename in matrix_files:
        assert os.path.exists(filename), "File not found: %s" % filename

    seqfile = None
    labeled_matrix_files = []
    outfiles = []
    try:
        # Write the sequences to a temp file.
        x, seqfile = tempfile.mkstemp(dir="."); os.close(x)
        handle = open(seqfile, 'w')
        for i, seq in enumerate(sequences):
            name = str(i)
            print >>handle, "%s \\%s\\" % (name, seq)
        handle.close()

        # Write the matrices to temp files.
        labeled_matrix_files = [_make_labeled_matrix(x) for x in matrix_files]

        # Set up temporary outfiles.
        for i in range(len(labeled_matrix_files)):
            x, f = tempfile.mkstemp(dir="."); os.close(x)
            outfiles.append(f)
        assert len(labeled_matrix_files) == len(outfiles)
        
        # Set up the patser commands to run.
        commands = []
        for i, mf in enumerate(labeled_matrix_files):
            x = _make_patser_cmd(patser_bin, seqfile, mf, alphabet_file)
            x = "%s >& %s" % (x, outfiles[i])
            commands.append(x)

        # Run patser either in serial or parallel.
        if num_jobs == 1 or len(commands) == 1 or not parallel.find():
            for cmd in commands:
                p = subprocess.Popen(
                    cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    close_fds=True)
                w, r = p.stdin, p.stdout
                w.close()
                r.read()
        else:
            parallel.run_commands(commands, num_jobs)

        # Parse the results.
        data = []
        for i, filename in enumerate(outfiles):
            #print open(filename).read(),; sys.exit(0)
            for x in parse(filename):
                name, position, strand, score, ln_p_value = x
                sequence_num = int(name)
                matrix_num = i
                # position is 1-based.  Convert to 0-based.
                position -= 1
                x = (sequence_num, matrix_num, position, strand, score,
                     ln_p_value)
                data.append(x)
    finally:
        files = [seqfile] + labeled_matrix_files + outfiles
        for x in files:
            if x and os.path.exists(x):
                os.unlink(x)

    return data

def patser(patser_bin, seq_file, matrix_file, alphabet_file, outhandle=None):
    import subprocess

    outhandle = outhandle or sys.stdout

    # Label the matrix_file first, and then pass to patser.
    file = None
    try:
        file = _make_labeled_matrix(matrix_file)
        cmd = _make_patser_cmd(patser_bin, seq_file, file, alphabet_file)
        p = subprocess.Popen(
            cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
        w, r = p.stdin, p.stdout
        w.close()
        for line in r:
            outhandle.write(line)
    finally:
        if file and os.path.exists(file):
            os.unlink(file)

def _make_labeled_matrix(matrix_file):
    # Return name of temporary file with labeled matrix.
    import time
    import tempfile
    
    # Name the file meaningfully, because it will show up in the
    # patser output.
    matrix_name = os.path.split(matrix_file)[-1]
    unique_str = ".%d.%g" % (os.getpid(), time.time())
    prefix = matrix_name + "."  # need "." to separate matrix from temp str
    suffix = unique_str
    x, file = tempfile.mkstemp(suffix=suffix, prefix=prefix, dir=".")
    os.close(x)
    #path = tempfile.gettempdir()
    #file = "%s/%s.%s" % (path, matrix_name, unique_str)

    x = label_matrix(matrix_file)
    open(file, 'w').write(x)

    return file

def _make_patser_cmd(patser_bin, seq_file, matrix_file, alphabet_file):
    cmd = [patser_bin]
    cmd.append("-f %s" % seq_file)
    cmd.append("-m %s" % matrix_file)
    cmd.append("-a %s" % alphabet_file)
    cmd.append("-c")     # search complementary strand
    cmd.append("-d2")    # ignore "N"
    cmd.append("-l 0")   # only print high scores (gone in patser-3e)
    if is_vertical_matrix(matrix_file):
        cmd.append("-v")
    if is_weight_matrix(matrix_file):
        cmd.append("-w")
    cmd = " ".join(cmd)
    return cmd

def parse(handle):
    # yields: name, position (1-based), strand, score, ln_p_value
    from filelib import openfh

    handle = openfh(handle)
    line = handle.readline()
    assert line.startswith("COMMAND LINE"), "Unexpected: %s" % line.strip()

    for line in handle:
        # Check for errors.
        if line.find("cannot execute binary file") >= 0:
            raise AssertionError, line.strip()
        # "-l" (item 9 on the command line) does not match any of the
        #    legal options.
        if line.find("does not match any of the legal options") >= 0:
            raise AssertionError, line.strip()
        
        # Bug: can discard useful information.
        if line.find("position") < 0:
            continue
        cols = line.strip().split()
        name, position, score, ln_p_value = \
              cols[0], cols[2], cols[4], cols[6]

        assert cols[1] == "position=" and cols[5] == "ln(p-value)="
        if position.endswith("C"):
            position = position[:-1]
            strand = "-"
        else:
            strand = "+"
        position = int(position)
        score, nlp = float(score), -float(ln_p_value)

        yield name, position, strand, score, nlp
        
def read_matrix(matrix_file):
    # Return horizontal matrix of counts or weights.  The rows
    # correspond to values for "A", "C", "G", "T", and the columns are
    # the positions in the matrix.
    vertical = is_vertical_matrix(matrix_file)
    labelled = is_labelled_matrix(matrix_file)
    
    lines = open(matrix_file).readlines()
    lines = [x.strip().split() for x in lines]

    if labelled:
        if vertical:
            lines = lines[1:]
        else:
            lines = [x[1:] for x in lines]

    # Convert to horizontal.
    if vertical:
        lines = [[x[i] for x in lines] for i in range(len(lines[0]))]
        
    fn = int
    if is_weight_matrix(matrix_file):
        fn = float
    lines = [[fn(x) for x in y] for y in lines]
    
    return lines

def label_matrix(matrix_file):
    # return a string
    if is_labelled_matrix(matrix_file):
        return open(matrix_file).read()
    
    lines = open(matrix_file).readlines()
    lines = [x.strip().split() for x in lines]

    # Assume the bases are ACGT.
    bases = ["A", "C", "G", "T"]
    if is_vertical_matrix(matrix_file):
        lines = [bases] + lines
    else:
        lines = [[bases[i]] + lines[i] for i in range(len(bases))]

    # Pretty format each column.
    widths = [0] * len(lines[0])
    for row in lines:
        for i in range(len(row)):
            widths[i] = max(widths[i], len(row[i]))
    outlines = []
    for row in lines:
        x = ["%*s" % (widths[i], row[i]) for i in range(len(row))]
        outlines.append(" ".join(x))
    return "\n".join(outlines) + "\n"

def is_vertical_matrix(matrix_file):
    # See if it's a vertical matrix.
    rows = open(matrix_file).readlines()
    num_rows = len(rows)
    
    num_cols = [len(x.strip().split()) for x in rows]
    assert min(num_cols) == max(num_cols), \
           "inconsistent columns in %s" % matrix_file
    num_cols = num_cols[0]

    # HACK: Is exactly 4 bases long, so ambiguous.
    if matrix_file.find("MA0094"):
        return 0
    assert not (num_cols == 4 and num_rows == 4), "ambiguous %s" % matrix_file
    
    if num_cols == 4:
        return 1
    elif num_rows == 4:
        return 0
    raise AssertionError, matrix_file

def is_weight_matrix(matrix_file):
    # HACK: Look for a decimal point anywhere.
    for x in open(matrix_file).readlines():
        if x.find(".") >= 0:
            return 1
    return 0

def is_labelled_matrix(matrix_file):
    lines = open(matrix_file).readlines()
    lines = [x.strip().split() for x in lines]

    if is_vertical_matrix(matrix_file):
        header = lines[0]
    else:
        header = [x[0] for x in lines]
    header.sort()
    header = "".join(header)
    return header == "ACGT"
