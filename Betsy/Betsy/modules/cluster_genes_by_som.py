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
        import cluster_genes_by_hierarchical

        filelib.safe_mkdir(out_path)
        metadata = {}

        raise NotImplementedError

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
        som_rows = mlib.get_user_option(
            user_options, "som_rows", not_empty=True, type=int)
        som_cols = mlib.get_user_option(
            user_options, "som_cols", not_empty=True, type=int)
        assert som_rows >= 1 and som_rows < 100
        assert som_cols >= 1 and som_cols < 100
        
        jobname = "cluster"
        cmd = cluster30.cluster30_file(
            in_data.identifier,
            (cluster_genes=="yes"), (cluster_arrays=="yes"), "som",
            distance=distance_metric, som_rows=som_rows, som_cols=som_cols,
            jobname=jobname)
        metadata["command"] = cmd

        # Find the output files and name them appropriately.
        cluster_files = cluster30._find_cluster_files(jobname)
        cluster_genes_by_hierarchical.fix_cluster30_dup_header(
            cluster_files["cdt"])

        opj = os.path.join
        out_cdt_file = opj(out_path, "signal.cdt")
        #out_kag_file = opj(out_path, "array_cluster.kag")
        #out_kgg_file = opj(out_path, "gene_cluster.kgg")

        assert "txt" in cluster_files
        shutil.copy2(cluster_files["txt"], out_cdt_file)
        #if "kag" in cluster_files:
        #    shutil.copy2(cluster_files["kag"], out_kag_file)
        #if "kgg" in cluster_files:
        #    shutil.copy2(cluster_files["kgg"], out_kgg_file)
        
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "cluster"
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #filename = 'cluster_file_' + original_file + '.cdt'
        #return filename
