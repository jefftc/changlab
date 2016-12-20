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
        from Betsy import module_utils as mlib
        import cluster_genes_by_hierarchical as clust
        
        filelib.safe_mkdir(out_path)
        metadata = {}

        kmeans_k = mlib.get_user_option(
            user_options, "kmeans_k", not_empty=True, type=int)
        assert kmeans_k >= 2 and kmeans_k < 100

        x = clust.run_cluster30(
            in_data.identifier, "kmeans", user_options, kmeans_k=kmeans_k)
        cmd, cluster_files = x
        metadata["command"] = cmd
        
        opj = os.path.join
        out_cdt_file = opj(out_path, "signal.cdt")
        out_kag_file = opj(out_path, "array_cluster.kag")
        out_kgg_file = opj(out_path, "gene_cluster.kgg")

        assert "cdt" in cluster_files
        shutil.copy2(cluster_files["cdt"], out_cdt_file)
        if "kag" in cluster_files:
            shutil.copy2(cluster_files["kag"], out_kag_file)
        if "kgg" in cluster_files:
            shutil.copy2(cluster_files["kgg"], out_kgg_file)
        
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "cluster"
