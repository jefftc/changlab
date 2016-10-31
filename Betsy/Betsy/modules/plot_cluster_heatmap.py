from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import filelib
        from Betsy import module_utils as mlib
        from plot_signal_heatmap import plot_heatmap

        metadata = {}
        cluster_files = mlib.find_cluster_files(in_data.identifier)
        assert "cdt" in cluster_files
        cmd = plot_heatmap(
            cluster_files["cdt"], outfile, cluster_files, user_options)
        metadata["command"] = cmd
        filelib.assert_exists_nz(outfile)
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "heatmap.png"
