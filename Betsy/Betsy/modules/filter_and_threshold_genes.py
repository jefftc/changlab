from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import filelib
        from genomicode import parallel
        from genomicode import config
        from genomicode import parselib
        from Betsy import module_utils as mlib
        from Betsy import bie3

        metadata = {}

        min_threshold = mlib.get_user_option(
            user_options, "min_threshold", type=float)
        assert not min_threshold or min_threshold >= -10
        max_threshold = mlib.get_user_option(
            user_options, "max_threshold", type=float)
        assert not max_threshold or max_threshold < 100000
        min_fold_change = mlib.get_user_option(
            user_options, "min_fold_change", type=float)
        assert not min_fold_change or (
            min_fold_change >= 0.1 and min_fold_change < 1000)
        min_delta = mlib.get_user_option(
            user_options, "min_delta", type=float)
        assert not min_delta or (min_delta >= 1 and min_delta < 100000)
        highest_var = mlib.get_user_option(
            user_options, "genes_with_highest_var", type=int)
        assert not highest_var or (highest_var >= 1 and highest_var < 20000)

        # Make a useful error message if no user options are
        # specified.
        # Find this module node.
        x = network.nodes
        x = [x for x in x if isinstance(x, bie3.ModuleNode)]
        x = [x for x in x if x.name == "filter_and_threshold_genes"]
        # This should not happen.
        assert x, "Missing module node"
        module_node = x[0]
        # Print out each of the options.
        msg = []
        for i, x in enumerate(module_node.option_defs):
            x = parselib.linesplit("%d.  %s.  %s" % (i+1, x.name, x.help))
            msg.extend(x)
        msg = "\n".join(msg)
        assert not (
            min_threshold is None and
            max_threshold is None and
            min_fold_change is None and
            min_delta is None and
            highest_var is None
            ), "Please specify criteria for filter_and_threshold:\n%s" % msg

        #if min_threshold:
        #    metadata["min_threshold"] = min_threshold
        #if max_threshold:
        #    metadata["max_threshold"] = max_threshold
        #if min_fold_change:
        #    metadata["min_fold_change"] = min_fold_change
        #if min_delta:
        #    metadata["min_delta"] = min_delta
        #if highest_var:
        #    metadata["highest_var"] = highest_var

        slice_matrix = filelib.which_assert(config.slice_matrix)

        sq = parallel.quote
        cmd = [
            sq(slice_matrix),
            ]
        if min_threshold:
            cmd += ["--min_value", min_threshold]
        if max_threshold:
            cmd += ["--max_value", max_threshold]
        if min_fold_change:
            #cmd += ["--select_row_fc_not_logged", min_fold_change]
            cmd += ["--select_row_fc", min_fold_change]
        if min_delta:
            cmd += ["--select_row_delta", min_delta]
        if highest_var:
            cmd += ["--select_row_var", highest_var]
        cmd += [sq(in_data.identifier)]
        cmd = " ".join(map(str, cmd))
        cmd = "%s > %s" % (cmd, sq(outfile))
        parallel.sshell(cmd)
        metadata["commands"] = [cmd]

        filelib.assert_exists_nz(outfile)
        return metadata
        
        #M = arrayio.read(in_data.identifier)
        #X = M.slice()
        #I_good = []
        #for i in range(M.nrow()):
        #    for j in range(len(X[i])):
        #        # Ignore missing values.
        #        if X[i][j] is None:
        #            continue
        #        if X[i][j] < min_threshold:
        #            M._X[i][j] = min_threshold
        #        if X[i][j] > max_threshold:
        #            M._X[i][j] = max_threshold
        #    x = M._X[i]
        #    x = [x for x in x if x is not None]
        #    if not x:
        #        # All missing values.  Ignore this gene.
        #        I_good.append(i)
        #        continue
        #    gene = x
        #    
        #    # what if min(gene) == 0?
        #    fold_change = max(gene) / float(min(gene))
        #    delta = max(gene) - min(gene)
        #    if fold_change >= min_fold_change and delta >= min_delta:
        #        I_good.append(i)
        # 
        #f = file(outfile, 'w')
        #M_c = M.matrix(I_good, None)
        #arrayio.tab_delimited_format.write(M_c, f)
        #f.close()

    def name_outfile(self, antecedents, user_options):
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #filename = 'signal_preprocessdataset_' + original_file + '.tdf'
        #return filename
        return "dataset.tdf"



