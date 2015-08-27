from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import arrayio
        from genomicode import arraysetlib

        signal_node1, signal_node2 = antecedents
        M1 = arrayio.read(signal_node1.identifier)
        M2 = arrayio.read(signal_node2.identifier)
        samples1 = M1.col_names(arrayio.COL_ID)
        samples2 = M2.col_names(arrayio.COL_ID)

        # Make sure no duplicate sample names.
        samples = samples1 + samples2
        seen = {}
        dups = {}
        for x in samples:
            if x in seen:
                dups[x] = 1
            seen[x] = 1
        dups = sorted(dups)
        assert not dups, "Duplicate sample names: %s" % ", ".join(dups)

        assert signal_node1.data.attributes["contents"] == "class0"
        assert signal_node2.data.attributes["contents"] == "class1"
        x1 = [0] * len(samples1)
        x2 = [1] * len(samples2)
        classes = x1 + x2
        arraysetlib.write_cls_file(outfile, "class0", "class1", classes)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        data_node1, data_node2 = antecedents
        original_file = module_utils.get_inputid(data_node1.identifier)
        filename = 'merge_' + original_file + '.cls'
        return filename


    
