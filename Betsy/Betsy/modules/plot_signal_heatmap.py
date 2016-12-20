from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import filelib

        metadata = {}
        cmd = plot_heatmap(in_data.identifier, outfile, {}, user_options)
        metadata["command"] = cmd

        #M = arrayio.read(in_data.identifier)
        #nrow = M.nrow()
        #ncol = M.ncol()
        #ratio = float(nrow) / ncol
        #max_box_height = 20
        #max_box_width = 60
    
        #if 'hm_width' in user_options:
        #    max_box_width = user_options['hm_width']
        #if 'hm_height' in user_options:
        #    max_box_height = user_options['hm_height']
        
        #if ratio >= 4:
        #    x, y = graphlib.find_tall_heatmap_size(
        #        nrow, ncol,
        #        max_box_height=max_box_height,
        #        max_box_width=max_box_width,
        #        min_box_height=20,
        #        min_box_width=20,
        #        max_megapixels=128)
        #else:
        #    x, y = graphlib.find_wide_heatmap_size(
        #        nrow, ncol,
        #        max_box_height=max_box_height,
        #        max_box_width=max_box_width,
        #        min_box_height=20,
        #        min_box_width=20,
        #        max_megapixels=128)
        #command.extend(['-x', str(x), '-y', str(y)])
        #
        #process = subprocess.Popen(command,
        #                           shell=False,
        #                           stdout=subprocess.PIPE,
        #                           stderr=subprocess.PIPE)
        #error_message = process.communicate()[1]
        #if error_message:
        #    raise ValueError(error_message)

        filelib.assert_exists_nz(outfile)
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "heatmap.png"


def plot_heatmap(filename, outfile, cluster_files, user_options):
    from genomicode import parallel
    from genomicode import filelib
    from genomicode import graphlib
    from Betsy import module_utils as mlib
    
    python = mlib.get_config(
        "python", which_assert_file=True, assert_exists=True)
    arrayplot = mlib.get_config(
        "arrayplot", which_assert_file=True, assert_exists=True)

    COLORS = [
        "red", "white", "red-green", "blue-yellow", "red-green-soft",
        "red-blue-soft", "matlab", "bild", "genepattern", "genespring",
        "yahoo", "brewer-prgn-div", "brewer-rdbu-div", 
        "brewer-rdylbu-div", "brewer-rdylgn-div", "brewer-spectral-div",
        "brewer-blues-seq", "brewer-greens-seq", "brewer-reds-seq",
        "brewer-ylorbr-seq", "brewer-qual-set",
        ]
    YESNO = ["no", "yes"]

    hm_width = mlib.get_user_option(user_options, "hm_width", type=int)
    hm_height = mlib.get_user_option(user_options, "hm_height", type=int)
    hm_color = mlib.get_user_option(
        user_options, "hm_color", allowed_values=COLORS, not_empty=True)

    hm_colorbar = mlib.get_user_option(
        user_options, "hm_colorbar", not_empty=True, allowed_values=YESNO)
    hm_colorbar_horizontal = mlib.get_user_option(
        user_options, "hm_colorbar_horizontal", not_empty=True,
        allowed_values=YESNO)
    hm_colorbar_height = mlib.get_user_option(
        user_options, "hm_colorbar_height", not_empty=True, type=float)
    hm_colorbar_width = mlib.get_user_option(
        user_options, "hm_colorbar_width", not_empty=True, type=float)
    hm_colorbar_font = mlib.get_user_option(
        user_options, "hm_colorbar_font", not_empty=True, type=float)

    hm_label_genes = mlib.get_user_option(
        user_options, "hm_label_genes", allowed_values=YESNO)
    hm_scale_gene_labels = mlib.get_user_option(
        user_options, "hm_scale_gene_labels", not_empty=True, type=float)
    hm_label_arrays = mlib.get_user_option(
        user_options, "hm_label_arrays", allowed_values=YESNO)
    hm_scale_array_labels = mlib.get_user_option(
        user_options, "hm_scale_array_labels", not_empty=True, type=float)

    hm_show_gene_tree = None
    hm_show_array_tree = None
    hm_show_gene_cluster = None
    hm_show_array_cluster = None
    if "hm_show_gene_tree" in user_options:
        hm_show_gene_tree = mlib.get_user_option(
            user_options, "hm_show_gene_tree", allowed_values=YESNO,
            not_empty=True)
        hm_show_array_tree = mlib.get_user_option(
            user_options, "hm_show_array_tree", allowed_values=YESNO,
            not_empty=True)
        hm_show_gene_cluster = mlib.get_user_option(
            user_options, "hm_show_gene_cluster", allowed_values=YESNO,
            not_empty=True)
        hm_show_array_cluster = mlib.get_user_option(
            user_options, "hm_show_array_cluster", allowed_values=YESNO,
            not_empty=True)

    # Set default values.
    if not hm_width or not hm_height:
        nrow, ncol = get_heatmap_size(filename)
        fn = graphlib.find_wide_heatmap_size
        if nrow > ncol:
            fn = graphlib.find_tall_heatmap_size
        x = fn(
            nrow, ncol, max_total_height=4096, max_total_width=4096,
            max_box_height=200, max_box_width=200)
        hm_width, hm_height = x

    if not hm_label_genes:
        nrow, ncol = get_heatmap_size(filename)
        hm_label_genes = "no"
        if nrow <= 50:
            hm_label_genes = "yes"
    if not hm_label_arrays:
        nrow, ncol = get_heatmap_size(filename)
        hm_label_arrays = "no"
        if ncol <= 50:
            hm_label_arrays = "yes"
    
        
    # Check values.
    assert hm_width >= 1 and hm_width <= 256, "Invalid width: %s" % hm_width
    assert hm_height >= 1 and hm_height <= 256, \
           "Invalid height: %s" % hm_height
    assert hm_scale_gene_labels > 0 and hm_scale_gene_labels < 10
    assert hm_scale_array_labels > 0 and hm_scale_array_labels < 10

    sq = parallel.quote
    cmd = [
        sq(python),
        sq(arrayplot),
        "--grid",
        "-x", hm_width,
        "-y", hm_height,
        "--color", hm_color,
        ]
    if hm_colorbar == "yes":
        cmd += [
            "--colorbar",
            "--cb_height", hm_colorbar_height,
            "--cb_width", hm_colorbar_width,
            "--cb_font", hm_colorbar_font,
            ]
        if hm_colorbar_horizontal == "yes":
            cmd += ["--cb_horizontal"]

    if hm_label_genes == "yes":
        cmd += [
            "--label_genes",
            "--scale_gene_labels", hm_scale_gene_labels,
            ]
    if hm_label_arrays == "yes":
        cmd += [
            "--label_arrays",
            "--scale_array_labels", hm_scale_array_labels,
            ]
    if hm_show_gene_tree == "yes" and "gtr" in cluster_files:
        cmd += ["--gene_tree_file", cluster_files["gtr"]]
    if hm_show_array_tree == "yes" and "atr" in cluster_files:
        cmd += ["--array_tree_file", cluster_files["atr"]]
    if hm_show_gene_cluster == "yes" and "kgg" in cluster_files:
        cmd += ["--gene_cluster_file", cluster_files["kgg"]]
    if hm_show_array_cluster == "yes" and "kag" in cluster_files:
        cmd += ["--array_cluster_file", cluster_files["kag"]]
    cmd += [
        sq(filename),
        sq(outfile),
        ]
    cmd = " ".join(map(str, cmd))
    parallel.sshell(cmd)

    return cmd


def _get_heatmap_size_h(filename):
    import arrayio
    MATRIX = arrayio.read(filename)
    nrow, ncol = MATRIX.nrow(), MATRIX.ncol()
    assert nrow and ncol
    return nrow, ncol

HM_SIZE_CACHE = {}
def get_heatmap_size(filename):
    global HM_SIZE_CACHE
    if filename not in HM_SIZE_CACHE:
        HM_SIZE_CACHE[filename] = _get_heatmap_size_h(filename)
    return HM_SIZE_CACHE[filename]
