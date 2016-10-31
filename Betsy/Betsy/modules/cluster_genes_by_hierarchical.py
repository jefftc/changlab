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

        DISTANCE_MEASURES = cluster30.DIST2ID.keys()
        LINKAGES = cluster30.METHOD2ID.keys()
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
        linkage = mlib.get_user_option(
            user_options, "linkage", not_empty=True,
            allowed_values=LINKAGES)

        jobname = "cluster"
        cmd = cluster30.cluster30_file(
            in_data.identifier,
            (cluster_genes=="yes"), (cluster_arrays=="yes"),
            "hierarchical", distance=distance_metric, method=linkage,
            jobname=jobname)
        metadata["command"] = cmd

        # Find the output files and name them appropriately.
        cluster_files = cluster30._find_cluster_files(jobname)

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
