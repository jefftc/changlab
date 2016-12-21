from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        data_node, rename_node = antecedents
        metadata = {}
        cmd = relabel(
            data_node.identifier, rename_node.identifier, outfile,
            user_options)
        metadata["command"] = cmd
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "signal.tdf"

    #def set_out_attributes(self, antecedents, out_attributes):
    #    data_node, rename_node = antecedents
    #    new_parameters = data_node.data.attributes.copy()
    #    new_parameters['relabel_sample'] = 'yes'
    #    return new_parameters


def relabel(data_file, rename_file, outfile, user_options):
    from genomicode import filelib
    from genomicode import parallel
    from Betsy import module_utils as mlib

    sample_header = mlib.get_user_option(
        user_options, "sample_labels_header", not_empty=True)
    # Make sure sample_header in rename file.
    x = open(rename_file).readline()
    x = x.rstrip("\r\n").split("\t")
    assert sample_header in x, "Missing header (%s): %s" % (
        sample_header, rename_file)

    sq = parallel.quote
    slice_matrix = mlib.get_config("slice_matrix", which_assert_file=True)
    x = "'%s,%s'" % (rename_file, sample_header)
    cmd = [
        "python",
        sq(slice_matrix),
        '--relabel_col_ids', x,
        sq(data_file),
        ]
    cmd = " ".join(cmd)
    cmd = "%s >& %s" % (cmd, outfile)
    parallel.sshell(cmd)

    filelib.assert_exists_nz(outfile)
    return cmd



