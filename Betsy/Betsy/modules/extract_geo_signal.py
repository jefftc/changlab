from Module import AbstractModule


class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import filelib

        outhandle = open(outfile, 'w')
        extract_signal(in_data.identifier, outhandle)
        outhandle.close()
        filelib.assert_exists_nz(outfile)

    def name_outfile(self, antecedents, user_options):
        return "signal.txt"



class FileMatrix:
    # BUG: Does not work on compressed files.
    def __init__(self, filename, coerce_fn=None):
        import os
        import re

        assert os.path.exists(filename)
        self.filename = filename
        self.coerce_fn = coerce_fn

        self.handle = None
        self.filemode = None

        # Collect some data about the size of the matrix.
        num_rows = num_cols = 0
        offset = 0
        offsets = []
        for line in self.get_handle():
            nc = len(re.findall("\t", line))+1
            if not num_cols:
                num_cols = nc
            assert num_cols == nc, "column mismatch: %s" % filename

            offsets.append(offset)
            num_rows += 1
            offset += len(line)
        self.num_rows, self.num_cols = num_rows, num_cols
        self.offsets = offsets

    def get_handle(self, mode="r"):
        if self.handle is None or self.filemode != mode:
            if self.handle is not None:
                self.handle.close()
            self.handle = open(self.filename, mode)
            self.filemode = mode
        return self.handle

    def close(self):
        if not self.handle:
            return
        self.handle.close()
        self.handle = None

    def __getitem__(self, index):
        handle = self.get_handle()
        assert self.num_rows == len(self.offsets)
        assert index < len(self.offsets), "%d %d" % (index, len(self.offsets))
        handle.seek(self.offsets[index])
        x = handle.readline()
        x = x.rstrip("\r\n").split("\t")
        if self.coerce_fn:
            x = [self.coerce_fn(x) for x in x]
        return x

    def __len__(self):
        assert self.num_rows == len(self.offsets)
        return self.num_rows

    def append(self, cols):
        handle = self.get_handle("a")
        handle.seek(0, 2)   # seek to end of file
        offset = handle.tell()
        print >>handle, "\t".join(map(str, cols))

        self.offsets.append(offset)
        self.num_rows += 1

    def __del__(self):
        self.close()

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]


def remove_quotes(s):
    if s and s[0] == '"' and s[-1] == '"':
        s = s[1:-1]
    return s

def extract_signal(filename, outhandle):
    import os
    import tempfile
    from genomicode import filelib

    # Write stuff to file to handle large data sets.
    tmpfile1 = tmpfile2 = tmpfile3 = None
    try:
        # tmpfile1        Raw signal data from series matrix file.
        # tmpfile2.<num>  Raw data split into separate tables.
        # tmpfile3        Final merged signal table.
        x, tmpfile1 = tempfile.mkstemp(dir=".")
        os.close(x)
        x, tmpfile2 = tempfile.mkstemp(dir=".")
        os.close(x)
        x, tmpfile3 = tempfile.mkstemp(dir=".")
        os.close(x)

        # Get a list of all lines in the series matrix tables.
        handle = open(tmpfile1, 'w')
        in_matrix_table = 0
        for cols in filelib.read_cols(filename):
            # Some files can have blank lines.
            if not cols:
                continue
            if cols[0] == "!series_matrix_table_begin":
                in_matrix_table = 1
            elif cols[0] == "!series_matrix_table_end":
                in_matrix_table = 0
            elif in_matrix_table:
                cols = [remove_quotes(x).strip() for x in cols]
                print >>handle, "\t".join(cols)
        handle.close()
        handle = None

        # Split the data into separate tables.
        num_tables = 0
        for line in filelib.openfh(tmpfile1):
            if line.startswith("ID_REF"):
                handle = open("%s.%d" % (tmpfile2, num_tables), 'w')
                num_tables += 1
            assert handle
            print >>handle, line,
        if handle:
            handle.close()
        assert num_tables

        # Sometimes the tables will not be aligned.
        # E.g. GSE9899-GPL570 contains two tables, and the 2nd is
        # missing some probe sets.  Get a list of the probe sets in
        # the tables.
        files = ["%s.%d" % (tmpfile2, i) for i in range(num_tables)]
        matrices = [FileMatrix(x) for x in files]
        id2indexes = []
        for matrix in matrices:
            id2index = {}
            for i, row in enumerate(matrix):
                id_ = row[0]
                id2index[id_] = i
            id2indexes.append(id2index)

        # Make a list of all the IDs.
        all_ids = {}
        for id2index in id2indexes:
            for id_ in id2index:
                all_ids[id_] = 1
        del all_ids["ID_REF"]
        all_ids = all_ids.keys()
        all_ids.sort()
        all_ids = ["ID_REF"] + all_ids

        # Align the indexes.
        #num_rows = row_names = None
        #for i in range(num_tables):
        #    filename = "%s.%d" % (tmpfile2, i)
        #    rname, nrow = [], 0
        #    for line in openfh(filename):
        #        x = line.split("\t", 1)[0]
        #        rname.append(x)
        #        nrow += 1
        #    if num_rows is None:
        #        num_rows = nrow
        #    if row_names is None:
        #        row_names = rname
        #    assert num_rows == nrow, "table is unaligned"
        #    assert row_names == rname

        # Merge all the pieces together into one big table.
        handle = open(tmpfile3, 'w')
        for id_ in all_ids:
            cols = []
            for matrix, id2index in zip(matrices, id2indexes):
                if id_ in id2index:
                    x = matrix[id2index[id_]]
                else:
                    # If this ID is missing, then just insert blank values.
                    x = [""] * len(matrix[0])
                if cols:
                    # If this is not the first matrix, then delete the
                    # row names.
                    x = x[1:]
                cols.extend(x)
            print >>handle, "\t".join(cols)
        handle.close()

        num_rows = len(all_ids)
        num_cols = len(filelib.read_cols(tmpfile3).next())

        # Figure out which expression values are missing.
        data_missing = {}
        for i, cols in enumerate(filelib.read_cols(tmpfile3)):
            assert len(cols) == num_cols, "line %d unaligned [%d:%d]" % (
                i, len(cols), num_cols)
            if i == 0:
                continue
            for j in range(1, len(cols)):
                try:
                    float(cols[j])
                except ValueError, x:
                    data_missing[(i, j)] = 1
                if cols[j] == "nan":
                    data_missing[(i, j)] = 1

        ## Remove the samples where >50% values are missing.
        #col_missing = [0] * num_cols   # number of values missing in each col
        #for i, j in data_missing:
        #    col_missing[j] += 1

        good_cols = [0]
        for i in range(1, num_cols):
            #if col_missing[i] > 0.50*(num_rows-1):  # -1 for the row names
            #    continue
            good_cols.append(i)

        ## Remove the genes where any value is missing.
        #row_missing = [0] * num_rows
        #for i, j in data_missing:
        #    if j not in good_cols:   # ignore samples that are already dropped
        #        continue
        #    row_missing[i] += 1

        good_rows = [0]
        for i in range(1, num_rows):
            #if row_missing[i] > 0:  # a value is missing.
            #    continue
            good_rows.append(i)

        assert len(good_cols) > 1, "no data"
        assert len(good_rows) > 1, "no data"

        # Write out the data.
        for i, cols in enumerate(filelib.read_cols(tmpfile3)):
            if i not in good_rows:
                continue
            x = [x for (i, x) in enumerate(cols) if i in good_cols]
            print >>outhandle, "\t".join(x)
    finally:
        x = [tmpfile1, tmpfile2, tmpfile3]
        x = [x for x in x if x]
        x = [x for x in x if os.path.exists(x)]
        for x in x:
            os.unlink(x)
        i = 0
        while 1:
            x = "%s.%d" % (tmpfile2, i)
            if os.path.exists(x):
                os.unlink(x)
            else:
                break
            i += 1
