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

        algorithm = out_attributes["de_algorithm"]
        metadata["algorithm"] = algorithm
        
        commands = calc_diffexp(
            data_node.identifier, cls_node.identifier, user_options, num_cores,
            out_path, algorithm)
        metadata["commands"] = commands
        metadata["num_cores"] = num_cores

        return metadata

    def name_outfile(self, antecedents, user_options):
        return "differential_expression"


def calc_diffexp(expression_file, class_label_file, user_options, num_cores,
                 out_path, algorithm):
    import os
    import arrayio
    from genomicode import parallel
    from genomicode import arraysetlib
    from genomicode import hashlib
    from genomicode import filelib
    from Betsy import module_utils as mlib

    fold_change_cutoff = mlib.get_user_option(
        user_options, "fold_change_cutoff", type=float)
    fdr_cutoff = mlib.get_user_option(
        user_options, "fdr_cutoff", type=float)
    p_cutoff = mlib.get_user_option(
        user_options, "p_cutoff", type=float)

    names, classes = arraysetlib.read_cls_file(class_label_file)
    assert names
    assert len(names) >= 2, (
        "At least 2 classes needed for differential expression "
        "analysis.  Found only: %s" % (names[0]))
    # Make sure there are the same number of samples in the class
    # label file as in the gene expression file.
    MATRIX = arrayio.read(expression_file)
    assert MATRIX.ncol() == len(classes), (
        "Mismatch: expression (%d) classes (%d)" % (
            MATRIX.ncol(), len(classes)))
    # Make sure classes go from [0, len(names))
    for i in classes:
        assert i >= 0 and i < len(names)

    x = []
    if fdr_cutoff:
        x.append("fdr_%g" % fdr_cutoff)
    if p_cutoff:
        x.append("p_%g" % p_cutoff)
    if fold_change_cutoff:
        x.append("fc_%g" % fold_change_cutoff)
    cutoff_str = ".".join(x)

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

            unfiltered_file = "%s.nocutoff.txt" % stem
            gmt_file = filtered_file = None
            if cutoff_str:
                gmt_file = "%s.%s.gmt" % (stem, cutoff_str)
                filtered_file = "%s.%s.txt" % (stem, cutoff_str)

            unfiltered_file = opj(out_path, unfiltered_file)
            if gmt_file:
                gmt_file = opj(out_path, gmt_file)
            if filtered_file:
                filtered_file = opj(out_path, filtered_file)

            # For heatmaps.
            clustered_file = "%s.%s.clustered.txt" % (stem, cutoff_str)
            row_tree_file = "%s.%s.row.txt" % (stem, cutoff_str)
            col_tree_file = "%s.%s.col.txt" % (stem, cutoff_str)
            heatmap_file = opj(
                out_path, "%s.%s.heatmap.png" % (stem, cutoff_str))
                
            x = filelib.GenericObject(
                N1=N1, N2=N2, I1=I1, I2=I2, stem=stem,
                fdr_cutoff=fdr_cutoff, p_cutoff=p_cutoff,
                fold_change_cutoff=fold_change_cutoff,
                unfiltered_file=unfiltered_file,
                filtered_file=filtered_file, gmt_file=gmt_file,
                clustered_file=clustered_file,
                row_tree_file=row_tree_file, col_tree_file=col_tree_file,
                heatmap_file=heatmap_file)
            jobs.append(x)

    commands = []
    for j in jobs:
        cmd = make_calc_diffexp_genes_command(
            expression_file, class_label_file,
            j.unfiltered_file, j.filtered_file, j.gmt_file,
            j.N1, j.N2, j.I1, j.I2, algorithm,
            j.fdr_cutoff, j.p_cutoff, j.fold_change_cutoff)
        commands.append(cmd)
    parallel.pshell(commands, max_procs=num_cores)

    # Make heatmaps for each of the filtered files.
    commands1 = []
    for j in jobs:
        cmd = make_cluster_command(
            j.filtered_file, j.clustered_file,
            j.row_tree_file, j.col_tree_file)
        commands1.append(cmd)
    commands1 = [x for x in commands1 if x]
    parallel.pshell(commands1, max_procs=num_cores)

    commands2 = []
    for j in jobs:
        cmd = make_heatmap_command(
            j.clustered_file, j.heatmap_file,
            j.row_tree_file, j.col_tree_file)
        commands2.append(cmd)
    commands2 = [x for x in commands2 if x]
    parallel.pshell(commands2, max_procs=num_cores)

    # Summarize the comparisons.
    if cutoff_str:
        fdr_cutoff_str = fdr_cutoff
        p_cutoff_str = p_cutoff
        fold_change_cutoff_str = fold_change_cutoff
        if fdr_cutoff_str is None:
            fdr_cutoff_str = ""
        if p_cutoff_str is None:
            p_cutoff_str = ""
        if fold_change_cutoff_str is None:
            fold_change_cutoff_str = ""
        summary_file = opj(out_path, "summary.%s.txt" % cutoff_str)
        handle = open(summary_file, 'w')
        header = [
            "Group 1", "Group 2", "FDR", "p", "Fold Change",
            "Higher in Group 1", "Higher in Group 2"]
        print >>handle, "\t".join(header)
        for j in jobs:
            filelib.assert_exists_nz(j.filtered_file)
            num_group_1 = num_group_2 = 0
            for d in filelib.read_row(j.filtered_file, header=1):
                HIGH = "Higher in "
                x = d.Direction
                assert x.startswith(HIGH)
                x = x[len(HIGH):].strip()
                assert x in [j.N1, j.N2]
                if x == j.N1:
                    num_group_1 += 1
                else:
                    num_group_2 += 1
            x = j.N1, j.N2, \
                fdr_cutoff_str, p_cutoff_str, fold_change_cutoff_str, \
                num_group_1, num_group_2
            assert len(x) == len(header)
            print >>handle, "\t".join(map(str, x))
        handle.close()
        mlib.txt2xls(summary_file, bold_header=True)
    
    return commands + commands1 + commands2


def make_calc_diffexp_genes_command(
    expression_file, class_label_file,
    unfiltered_file, filtered_file, gmt_file,
    name1, name2, indexes1, indexes2, algorithm,
    fdr_cutoff, p_cutoff, fold_change_cutoff):
    # indexes should be 1-based, not including headers.
    from genomicode import parallel
    from genomicode import filelib
    from genomicode import parselib
    from Betsy import module_utils as mlib
    from Betsy.rules import DiffExp

    filelib.assert_exists_nz(expression_file)
    filelib.assert_exists_nz(class_label_file)
    assert algorithm in DiffExp.DE_ALGORITHMS

    diffexp = mlib.get_config("find_diffexp_genes", which_assert_file=True)

    ranges1 = [(i, i+1) for i in indexes1]
    ranges2 = [(i, i+1) for i in indexes2]
    indexes1_str = parselib.unparse_ranges(ranges1)
    indexes2_str = parselib.unparse_ranges(ranges2)

    sq = parallel.quote
    cmd = [
        'python',
        sq(diffexp),
        sq(expression_file),
        '--algorithm', algorithm,
        "--name1", name1,
        "--name2", name2,
        "--indexes1", indexes1_str,
        "--indexes2", indexes2_str,
        ]
    if fdr_cutoff is not None:
        cmd += ["--fdr_cutoff", fdr_cutoff]
    if p_cutoff is not None:
        cmd += ["--p_cutoff", p_cutoff]
    if fold_change_cutoff is not None:
        cmd += ["--fold_change", fold_change_cutoff]
    # By default, save the unfiltered output.
    outfile = unfiltered_file
    if gmt_file and filtered_file:
        # If want filtering, move the files around.
        cmd += ["--gmt_file", sq(gmt_file)]
        cmd += ["--unfiltered_file", sq(unfiltered_file)]
        outfile = filtered_file

    cmd = " ".join(map(str, cmd))
    cmd = "%s >& %s" % (cmd, outfile)
    return cmd


def make_cluster_command(in_file, out_file, row_tree_file, col_tree_file):
    import arrayio
    from genomicode import filelib
    from genomicode import parallel
    from Betsy import module_utils as mlib

    # If this file is too big, then don't bother clustering.
    MAX_ROWS = 2500
    
    MATRIX = arrayio.read(in_file)
    if MATRIX.nrow() > MAX_ROWS:
        return None

    filelib.assert_exists_nz(in_file)
    
    slice_matrix = mlib.get_config("slice_matrix", which_assert_file=True)
    sq = parallel.quote
    cmd = [
        sq(slice_matrix),
        "--gc", "mean",
        "--gn", "var",
        "--reorder_col_cluster",
        "--col_tree_file", sq(col_tree_file),
        "--reorder_row_cluster",
        "--row_tree_file", sq(row_tree_file),
        sq(in_file),
        ]
    cmd = "%s > %s" % (" ".join(cmd), sq(out_file))
    return cmd


def make_heatmap_command(in_file, heatmap_file, row_tree_file, col_tree_file):
    import os
    from genomicode import parallel
    from Betsy import module_utils as mlib
    import plot_signal_heatmap

    if not os.path.exists(in_file):
        return None

    # Figure out dimensions of heatmap.
    width, height = plot_signal_heatmap.get_heatmap_size(in_file)

    arrayplot = mlib.get_config("arrayplot", which_assert_file=True)
    sq = parallel.quote
    cmd = [
        sq(arrayplot),
        "-x", width,
        "-y", height,
        "--grid",
        "--array_tree_file", sq(col_tree_file),
        "--al",
        "--gene_tree_file", sq(row_tree_file),
        "--gl",
        #"--color", "brewer-greens-seq",
        sq(in_file),
        sq(heatmap_file),
        ]
    cmd = " ".join(map(str, cmd))
    return cmd
    

    
    
    
