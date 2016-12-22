from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os

        data_node, cls_node = antecedents
        if not os.path.exists(out_path):
            os.mkdir(out_path)
        metadata = {}

        permutation_type = out_attributes["permutation_type"]
        database = out_attributes["geneset_database"]

        x = calc_gsea(
            data_node.identifier, cls_node.identifier, user_options, num_cores,
            out_path, permutation_type, database)
        commands, permutation_types = x
        permutation_types = ",".join(permutation_types)
        metadata["commands"] = commands
        metadata["permutation_type"] = permutation_types
        metadata["num_cores"] = 1

        return metadata


    def name_outfile(self, antecedents, user_options):
        return "gsea"


    
def calc_gsea(
    expression_file, class_label_file, user_options, num_cores, out_path,
    permutation_type, database):
    import os
    import arrayio
    from genomicode import parallel
    from genomicode import arraysetlib
    from genomicode import hashlib
    from genomicode import filelib

    names, classes = arraysetlib.read_cls_file(class_label_file)
    assert names
    assert len(names) >= 2, (
        "At least 2 classes needed for GSEA analysis.  "
        "Found only: %s" % (names[0]))
    # Make sure there are the same number of samples in the class
    # label file as in the gene expression file.
    MATRIX = arrayio.read(expression_file)
    assert MATRIX.ncol() == len(classes), (
        "Mismatch: expression (%d) classes (%d)" % (
            MATRIX.ncol(), len(classes)))
    # Make sure classes go from [0, len(names))
    for i in classes:
        assert i >= 0 and i < len(names)

    # Find all combinations of names and classes.
    opj = os.path.join
    jobs = []
    for i1 in range(len(names)-1):
        for i2 in range(i1+1, len(names)):
            N1 = names[i1]
            N2 = names[i2]
            # Indexes should be 1-based.
            I1 = [i+1 for i in range(len(classes)) if classes[i] == i1]
            I2 = [i+1 for i in range(len(classes)) if classes[i] == i2]
            N1_h = hashlib.hash_var(N1)
            N2_h = hashlib.hash_var(N2)
            stem = "%s.vs.%s" % (N1_h, N2_h)

            gsea_path = opj(out_path, "%s.%s.gsea" % (stem, database))
                
            x = filelib.GenericObject(
                N1=N1, N2=N2, I1=I1, I2=I2, stem=stem,
                gsea_path=gsea_path)
            jobs.append(x)

    permutation_types = {}
    commands = []
    for j in jobs:
        # Need at least 3 samples for "phenotype" permutations.  If
        # there are fewer samples, then set to "gene_set".
        if len(I1) < 3 or len(I2) < 3:
            permutation_type = "gene_set"
        permutation_types[permutation_type] = 1
        cmd = make_gsea_command(
            expression_file, class_label_file, j.gsea_path,
            j.N1, j.N2, j.I1, j.I2, permutation_type, database)
        commands.append(cmd)
    for cmd in commands:
        parallel.sshell(cmd)
    return commands, sorted(permutation_types)
    

def make_gsea_command(
    expression_file, class_label_file, gsea_path,
    name1, name2, indexes1, indexes2, permutation_type, database):
    # indexes should be 1-based, not including headers.
    from genomicode import parallel
    from genomicode import filelib
    from genomicode import parselib
    from Betsy import module_utils as mlib
    from Betsy.rules import GSEAAnalysis

    filelib.assert_exists_nz(expression_file)
    filelib.assert_exists_nz(class_label_file)
    assert permutation_type in GSEAAnalysis.GSEA_PERMUTATION
    assert database in GSEAAnalysis.GSEA_DATABASE

    ranges1 = [(i, i+1) for i in indexes1]
    ranges2 = [(i, i+1) for i in indexes2]
    indexes1_str = parselib.unparse_ranges(ranges1)
    indexes2_str = parselib.unparse_ranges(ranges2)

    gsea = mlib.get_config("gsea", which_assert_file=True)
    
    sq = parallel.quote
    cmd = [
        sq(gsea),
        "--name1", name1,
        "--name2", name2,
        "--indexes1", indexes1_str,
        "--indexes2", indexes2_str,
        "--permutation_type", sq(permutation_type),
        "--database", sq(database),
        "--min_match_score", 0.80,
        sq(expression_file),
        sq(gsea_path),
        ]
    cmd = " ".join(map(str, cmd))
    return cmd
