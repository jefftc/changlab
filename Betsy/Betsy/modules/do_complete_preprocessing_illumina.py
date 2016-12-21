from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        import shutil
        from genomicode import filelib

        #(data_node1, data_node2, data_node3, data_node4, data_node5,
        # data_node6, data_node7, data_node8, data_node9) = antecedents
        data_nodes = [
            ("SignalFile", "gene_expression.nonorm.txt"),
            ("SignalFile", "gene_expression.normalized.txt"),
            ("SignalDistributionBoxplot", "signal_distribution.png"),
            ("ActbPlot", "ACTB.nonorm.png"),       # beta-actin expression
            ("ActbPlot", "ACTB.normalized.png"),   # beta-actin expression
            ("PCAPlot", "pca.nonorm.png"),         # No normalization
            ("PCAPlot", "pca.normalized.png"),     # Normalized
            ("Heatmap", "heatmap.nonorm.png"),
            ("Heatmap", "heatmap.normalized.png"),
            ("IlluminaControlFile", "control_probes.gct"),  # Control probes.
            ("HousekeepingPlot", "housekeeping.expression.png"), 
            ("BiotinPlot", "biotin.control.png"),           # Biotin controls.
            ("IlluHybridizationProbePlot", "hybridization.control.png"),
            ]
        assert len(antecedents) == len(data_nodes)
        for i, (dtype, outfile) in enumerate(data_nodes):
            inode = antecedents[i]
            filelib.assert_exists_nz(inode.identifier)
            assert inode.data.datatype.name == dtype, "Mismatch: %s %s" % (
                inode.data.datatype.name, dtype)
        if not os.path.exists(out_path):
            os.mkdir(out_path)
        metadata = {}

        # Copy the files over.
        for i, (dtype, outfile) in enumerate(data_nodes):
            inode = antecedents[i]
            outfilename = os.path.join(out_path, outfile)
            shutil.copy2(inode.identifier, outfilename)

        return metadata

    
    def name_outfile(self, antecedents, user_options):
        return "report"


    
