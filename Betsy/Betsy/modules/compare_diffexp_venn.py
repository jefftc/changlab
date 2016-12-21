from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import genesetlib
        from genomicode import parallel
        from genomicode import filelib
        from Betsy import module_utils as mlib
        
        in_data = antecedents
        if not os.path.exists(out_path):
            os.mkdir(out_path)
        metadata = {}

        merge = mlib.get_user_option(
            user_options, "merge_up_and_down_genes",  not_empty=True,
            allowed_values=["yes", "no"])
        merge = (merge == "yes")

        opj = os.path.join
        gs_filename = opj(out_path, "gene_sets.gmt")
        intersect_filename = opj(out_path, "intersection.gmt")
        count_filename = opj(out_path, "pairwise_count_matrix.txt")
        venn_plot_file = opj(out_path, "venn.tiff")

        # Make a list of all the data sets in the antecedents.
        # <stem>.nocutoff.txt
        # <stem>.<cutoff>.txt
        # <stem>.<cutoff>.gmt
        # <stem>.<cutoff>.heatmap.png
        x = os.listdir(in_data.identifier)
        x = sorted(x)
        x = [x for x in x if x.endswith(".gmt")]
        x = [x for x in x if x.find("nocuttof") < 0]
        x = [opj(in_data.identifier, x) for x in x]
        filtered_geneset_files = x
        assert filtered_geneset_files, "Missing: filtered geneset files"

        # For each of the filtered_geneset_files, figure out the
        # <stem>.  This is tricky because <cutoff> may contain
        # multiple dots.
        # <stem>.fdr_0.05.p_0.05.fc_1.5.gmt
        stems = [None] * len(filtered_geneset_files)
        for i, x in enumerate(filtered_geneset_files):
            x = filtered_geneset_files[i]
            x = os.path.split(x)[1]
            x = x.split(".")
            j = 0
            while j < len(x):
                if x[j].startswith("fdr_") or x[j].startswith("fc_") or \
                       x[j].startswith("p_") or x[j] == "gmt":
                    x = x[:j]
                else:
                    j += 1
            x = ".".join(x)
            stems[i] = x

        genesets = []
        geneset_stems = []
        for i, filename in enumerate(filtered_geneset_files):
            for x in genesetlib.read_genesets(filename):
                name, description, genes = x
                x = genesetlib.GeneSet(name, description, genes)
                genesets.append(x)
                geneset_stems.append(stems[i])
        assert genesets, "I could not find any gene sets"
                
        # Should contain gene sets whose name fits the pattern:
        # <NAME>_ID_UP
        # <NAME>_ID_DN
        # <NAME>_NAME_UP
        # <NAME>_NAME_DN
        # Only want the _ID_ gene sets for comparison.
        I = [i for (i, x) in enumerate(genesets) if x.name.find("_ID_") >= 0]
        genesets = [genesets[i] for i in I]
        geneset_stems = [geneset_stems[i] for i in I]
        assert genesets, "I could not find any '_ID_' gene sets"

        # Rename each of the gene sets.
        for i, gs in enumerate(genesets):
            name = gs.name
            n = "%s_%s" % (geneset_stems[i], name)
            if merge:
                # If I'm merging, then which gene set is UP or DN
                # doesn't matter.
                suffix = name[-6:]
                assert suffix in ["_ID_UP", "_ID_DN"]
                n = "%s%s" % (geneset_stems[i], suffix)
            gs.name = n

        # Write out the gene sets.
        genesetlib.write_gmt(gs_filename, genesets)

        # Count the number of gene sets.
        x = [x.name for x in genesets]
        if merge:
            x = [x.replace("_DN", "_UP") for x in x]
            x = [x.replace("_DOWN", "_UP") for x in x]
        num_genesets = len({}.fromkeys(x))

        calc_venn = mlib.get_config("calc_venn", which_assert_file=True)
        sq = parallel.quote
        cmd = [
            sq(calc_venn),
            "-o", sq(intersect_filename),
            "--all_genesets",
            "--num_to_compare", 2,
            ]
        if num_genesets <= 5:
            # Can only plot up to 5 circles.
            cmd += ["--plotfile", sq(venn_plot_file)]
        if merge:
            cmd += ["--automatch"]
        cmd.append(sq(gs_filename))
        cmd = " ".join(map(str, cmd))
        cmd = "%s >& %s" % (cmd, sq(count_filename))
        parallel.sshell(cmd)
        metadata["commands"] = [cmd]

        # Make a heatmap of the counts.
        UNCLUSTERED_FILE = "unclustered.txt"
        CLUSTERED_FILE = "clustered.txt"
        COL_TREE_FILE = "col_tree.txt"
        ROW_TREE_FILE = "row_tree.txt"
        HEATMAP_FILE = opj(out_path, "heatmap.counts.png")

        # Make a file with the counts.
        outhandle = open(UNCLUSTERED_FILE, 'w')
        for line in open(count_filename):
            if not line.strip():
                break
            outhandle.write(line)
        outhandle.close()

        # Cluster the counts.
        slice_matrix = mlib.get_config("slice_matrix", which_assert_file=True)
        arrayplot = mlib.get_config("arrayplot", which_assert_file=True)
        cmd = [
            sq(slice_matrix),
            "--reorder_col_cluster",
            "--col_tree_file", sq(COL_TREE_FILE),
            "--reorder_row_cluster",
            "--row_tree_file", sq(ROW_TREE_FILE),
            sq(UNCLUSTERED_FILE),
            ]
        cmd = "%s > %s" % (" ".join(cmd), sq(CLUSTERED_FILE))
        parallel.sshell(cmd)
        metadata["commands"].append(cmd)
        filelib.assert_exists_nz(CLUSTERED_FILE)

        # Draw the heatmap.
        cmd = [
            sq(arrayplot),
            "--grid",
            "--array_tree_file", sq(COL_TREE_FILE),
            "--al",
            "--gene_tree_file", sq(ROW_TREE_FILE),
            "--gl",
            "--colorbar",
            "--color", "brewer-greens-seq",
            sq(CLUSTERED_FILE),
            sq(HEATMAP_FILE),
            ]
        cmd = " ".join(cmd)
        parallel.sshell(cmd)
        metadata["commands"].append(cmd)
        filelib.assert_exists_nz(HEATMAP_FILE)

        mlib.txt2xls(count_filename)

        return metadata

    def name_outfile(self, antecedents, user_options):
        return "venn"
