from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        import os
        from genomicode import filelib

        count_path = in_data.identifier
        assert os.path.exists(count_path)
        assert os.path.isdir(count_path)
        result_files = filelib.list_files_in_path(
            count_path, endswith=".count")
        assert result_files, "No .count files found."

        # Parse the count files.
        name2results = {}
        for filename in result_files:
            x = os.path.split(filename)[1]
            x = os.path.splitext(x)[0]
            name = x
            assert name not in name2results
            x = _read_htseq_count_file(filename)
            name2results[name] = x
        assert name2results, "No samples"

        # Assemble into a summary matrix.
        # Rows:
        # no_feature
        # ambiguous
        # too_low_aQual
        # not_aligned
        # alignment_not_unique
        # total_mapped
        # total_fragments
        # percent_mapped

        ROWS = [
            "no_feature",
            "ambiguous",
            "too_low_aQual",
            "not_aligned",
            "alignment_not_unique",
            ]
        all_names = sorted(name2results)
        
        matrix = []
        header = ["Feature"] + all_names
        matrix.append(header)
        for rn in ROWS:
            x = [rn] + [getattr(name2results[n], rn) for n in all_names]
            assert len(x) == len(header)
            matrix.append(x)
            
        # Count the total mapped and total_fragments.
        total_mapped = []
        total_fragments = []
        perc_mapped = []
        for n in all_names:
            # Sum up the counts
            results = name2results[n] 
            x1 = sum(results.counts.values())
            x2 = 0
            for rn in ROWS:
                x2 += getattr(results, rn)
            total_mapped.append(x1)
            total_fragments.append(x1+x2)
            perc_mapped.append(x1/float(x1+x2))

        x1 = ["total_mapped"] + total_mapped
        x2 = ["total_fragments"] + total_fragments
        x3 = ["perc_mapped"] + perc_mapped
        assert len(x1) == len(header)
        assert len(x2) == len(header)
        assert len(x3) == len(header)
        matrix.append(map(str, x1))
        matrix.append(map(str, x2))
        matrix.append(map(str, x3))

        
        # Write out the data file.
        handle = open(outfile, 'w')
        for x in matrix:
            print >>handle, "\t".join(map(str, x))


    def name_outfile(self, antecedents, user_options):
        return "summary.txt"


# Move into its own file.
class HTSeqCountResults:
    def __init__(
        self, counts, no_feature, ambiguous, too_low_aQual,
        not_aligned, alignment_not_unique, warnings, errors):
        # counts is a dictionary of gene_id -> count.
        # warnings is a list of warning lines
        # errors is a list of error lines
        self.counts = counts.copy()
        self.no_feature = no_feature        # count not be assigned to feature
        self.ambiguous = ambiguous          # multiple features
        self.too_low_aQual = too_low_aQual  # low quality
        self.not_aligned = not_aligned      # no alignment in SAM file
        self.alignment_not_unique = alignment_not_unique  # multiple alignment
        self.warnings = warnings[:]
        self.errors = errors[:]


def _read_htseq_count_file(file_or_handle):
    # Return an HTSeqCountResults object.
    import os

    #filename = None
    handle = file_or_handle
    if type(file_or_handle) is type(""):
        #filename = file_or_handle
        assert os.path.exists(file_or_handle)
        handle = open(file_or_handle)

    counts = {}
    meta = {}
    warnings = []
    errors = []
    for line in handle:
        # 100000 GFF lines processed.
        if line.rstrip().endswith("processed."):
            continue
        # Warning: Mate records missing for 128 records; first such
        # record: <SAM_Alignment object: Paired-end read
        # 'SRR988443.108873869' aligned to MT:[11627,11677)/->.
        if line.startswith("Warning:"):
            warnings.append(line.strip())
            continue
        #err_msg = line.strip()
        #if filename:
        #    err_msg = "%s: %s" % line.strip()
        #assert not line.startswith("Error occured"), err_msg
        if line.startswith("Error"):
            errors.append(line.strip())
            continue
        x = line.rstrip("\r\n").split("\t")
        assert len(x) == 2, "Unrecognized line: %s" % repr(line.rstrip("\r\n"))
        gene_id, count = x
        gene_id, count = gene_id.strip(), int(count)
        # __no_feature    343859
        # __ambiguous     279435
        # __too_low_aQual 247071
        # __not_aligned   9744
        # __alignment_not_unique  0
        if gene_id.startswith("__"):
            name = gene_id[2:]
            meta[name] = count
            continue
        counts[gene_id] = count
    assert len(meta) == 5

    x = HTSeqCountResults(counts, warnings=warnings, errors=errors, **meta)
    return x
