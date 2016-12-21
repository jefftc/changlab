from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        import shutil
        from genomicode import filelib
        from genomicode import cluster30
        from Betsy import module_utils as mlib

        filelib.safe_mkdir(out_path)
        metadata = {}
        
        LINKAGES = cluster30.METHOD2ID.keys()
        linkage = mlib.get_user_option(
            user_options, "linkage", not_empty=True,
            allowed_values=LINKAGES)

        x = run_cluster30(
            in_data.identifier, "hierarchical", user_options, method=linkage)
        cmd, cluster_files = x
        metadata["command"] = cmd

        opj = os.path.join
        out_cdt_file = opj(out_path, "signal.cdt")
        out_atr_file = opj(out_path, "array_tree.atr")
        out_gtr_file = opj(out_path, "gene_tree.gtr")

        assert "cdt" in cluster_files
        shutil.copy2(cluster_files["cdt"], out_cdt_file)
        if "atr" in cluster_files:
            shutil.copy2(cluster_files["atr"], out_atr_file)
        if "gtr" in cluster_files:
            shutil.copy2(cluster_files["gtr"], out_gtr_file)
        
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "cluster"
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #filename = "cluster_file_" + original_file + ".cdt"
        #return filename


def fix_cluster30_dup_header(filename):
    # Cluster30 creates a file with "NAME" as the header for the third
    # column.  If the infile also has a "NAME" column, then this will
    # be duplicated.  Detect this situation and fix it.
    from genomicode import filelib
    
    filelib.assert_exists_nz(filename)
    matrix = [x for x in filelib.read_cols(filename)]
    assert matrix
    assert matrix[0]
    header = matrix[0]
    # GID  <COL0>  NAME  GWEIGHT  <COL1>  [<SAMPLES>...]
    assert len(header) >= 5
    changed = False
    if header[1] == "NAME" and header[2] == "NAME":
        header[1] = "NAME_"
        changed = True
    if not changed:
        return
    handle = open(filename, 'w')
    for x in matrix:
        print >>handle, "\t".join(x)
    
    
def run_cluster30(filename, algorithm, user_options, **more_args):
    import arrayio
    from genomicode import cluster30
    from Betsy import module_utils as mlib

    MATRIX_FILE = "data.pcl"

    DISTANCE_MEASURES = cluster30.DIST2ID.keys()
    YESNO = ["yes", "no"]

    cluster_genes = mlib.get_user_option(
        user_options, "cluster_genes", not_empty=True,
        allowed_values=YESNO)
    cluster_arrays = mlib.get_user_option(
        user_options, "cluster_arrays", not_empty=True,
        allowed_values=YESNO)
    distance_metric = mlib.get_user_option(
        user_options, "distance_measure", not_empty=True,
        allowed_values=DISTANCE_MEASURES)

    # Make a PCL-formatted file for cluster 3.0.  It might
    # misinterpret the columns of a tab-delimited file.
    matrix = arrayio.read(filename)
    matrix = arrayio.convert(matrix, to_format=arrayio.pcl_format)
    arrayio.write(matrix, open(MATRIX_FILE, 'w'))

    jobname = "cluster"
    cmd = cluster30.cluster30_file(
        MATRIX_FILE, 
        (cluster_genes=="yes"), (cluster_arrays=="yes"),
        algorithm, distance=distance_metric, 
        jobname=jobname, **more_args)

    # Find the output files and name them appropriately.
    cluster_files = cluster30._find_cluster_files(jobname)
    fix_cluster30_dup_header(cluster_files["cdt"])

    return cmd, cluster_files
