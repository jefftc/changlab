from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_filename):
        import os
        from genomicode import filelib
        from genomicode import alignlib

        count_path = in_data.identifier
        assert os.path.exists(count_path)
        assert os.path.isdir(count_path)

        result_files = filelib.list_files_in_path(
            count_path, endswith=".count")
        assert result_files, "No .count files found."

        preprocess = out_attributes.get("preprocess")
        assert preprocess in ["counts"]

        # For each of the count files, read the count data.
        name2results = {}
        for filename in result_files:
            x = os.path.split(filename)[1]
            x = os.path.splitext(x)[0]
            name = x
            assert name not in name2results
            
            x = alignlib.parse_htseq_count_output(filename)
            name2results[name] = x
        assert name2results, "No samples"

        # Check for errors.
        has_errors = [x for x in name2results if name2results[x].errors]
        x = ""
        if len(has_errors) > 1:
            x = "s"
        assert not has_errors, "Errors found in sample%s: %s" % (
            x, ", ".join(has_errors))

        # Assemble into a gene expression matrix.
        all_names = sorted(name2results)
        all_genes = {}
        for n in all_names:
            all_genes.update(name2results[n].counts)
        all_genes = sorted(all_genes)
        # Get rid of <blank>.
        all_genes = [x for x in all_genes if x]

        header = ["Gene ID"] + all_names
        matrix = []
        matrix.append(header)
        for gene_id in all_genes:
            counts = [
                name2results[n].counts.get(gene_id, 0) for n in all_names]
            x = [gene_id] + counts
            matrix.append(x)
        
        # Write out the data file.
        handle = open(out_filename, 'w')
        for x in matrix:
            print >>handle, "\t".join(map(str, x))


    def name_outfile(self, antecedents, user_options):
        return "data.counts"


## def _read_htseq_count_file(file_or_handle):
##     # Return dictionary of <gene_id> -> <count>.
##     import os

##     filename = None
##     handle = file_or_handle
##     if type(file_or_handle) is type(""):
##         filename = file_or_handle
##         assert os.path.exists(file_or_handle)
##         handle = open(file_or_handle)

##     counts = {}
##     for line in handle:
##         # 100000 GFF lines processed.
##         if line.rstrip().endswith("processed."):
##             continue
##         # Warning: Mate records missing for 128 records; first such
##         # record: <SAM_Alignment object: Paired-end read
##         # 'SRR988443.108873869' aligned to MT:[11627,11677)/->.
##         if line.startswith("Warning: Mate records missing"):
##             continue
##         err_msg = line.strip()
##         if filename:
##             err_msg = "ERROR: %s" % line.strip()
##         assert not line.startswith("Error occured"), err_msg
##         # __no_feature    343859
##         # __ambiguous     279435
##         # __too_low_aQual 247071
##         # __not_aligned   9744
##         # __alignment_not_unique  0
##         if line.startswith("__"):
##             continue
##         x = line.rstrip("\r\n").split("\t")
##         assert len(x) == 2, "Unrecognized line: %s" % repr(line.rstrip("\r\n"))
##         gene_id, count = x
##         counts[gene_id.strip()] = int(count)
##     return counts
