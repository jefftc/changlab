from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import arrayio
        from genomicode import filelib
        from genomicode import arraysetlib
        #from Betsy import read_label_file
        #from Betsy import module_utils
        from Betsy.rules import BasicDataTypes as BDT

        signal_node, slf_node = antecedents
        slf_filename = slf_node.identifier

        samples = []  # Ordered list of sample names.
        classes = []  # Ordered list of classes
        for d in filelib.read_row(slf_filename, header=1):
            x = hasattr(d, "Sample") and hasattr(d, "Class")
            assert x, "SimpleLabelFile must have headers Sample and Class."
            samples.append(d.Sample)
            classes.append(d.Class)

        # Make into a dictionary.  Will lose information if sample
        # names not unique.
        sample2class = {}
        for i in range(len(samples)):
            assert samples[i] not in sample2class, \
                   "Duplicate sample: %s" % samples[i]
            sample2class[samples[i]] = classes[i]
        # Get a list of the unique classes.  Maintain order.
        class_names = []
        for x in classes:
            if x not in class_names:
                class_names.append(x)

        # Check whether number of classes matches contents variable.
        x = slf_node.data.attributes["contents"]
        num_classes = BDT.CONTENTS2NUMCLASSES[x]
        assert len(class_names) == num_classes, \
               'For files with "%s", I expected %d classes but found %d.' % (
            x, num_classes, len(class_names))
        
        # Get the names of the samples in the matrix.
        M = arrayio.read(signal_node.identifier)
        samples = M.col_names('_SAMPLE_NAME')
        # Make sure all my samples are known.
        for x in samples:
            assert x in sample2class, "Missing sample: %s" % x
        classes = [sample2class[x] for x in samples]

        arraysetlib.write_multi_cls_file(outfile, class_names, classes)


    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        signal_node, label_node = antecedents
        original_file = module_utils.get_inputid(signal_node.identifier)
        filename = 'class_label_' + original_file + '.cls'
        return filename
