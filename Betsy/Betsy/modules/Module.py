class AbstractModule:
    # Methods:
    # name_outfile        REQUIRED.
    # run                 REQUIRED.
    # set_out_attributes  OPTIONAL.  Needed if the module determines the
    #                     output attributes, e.g. is_logged.
    
    # DEPRECATE THIS.  MOVE AWAY
    # hash_input          OPTIONAL.  Maybe never need to overload?

    def __init__(self):
        pass


    def run(self, network, antecedents, out_attributes, user_options,
            num_cores, outfile):
        # Performs the computation that the module is supposed to do.
        # Does not return anything.  Can raise an exception of there
        # is a problem.
        #
        # There is no default implementation.
        #
        # network         Network object.  Likely not used.  But can be
        #                 useful when generating a report that shows the
        #                 structure of the network.
        # antecedents     Either a single or list of DataObjects, in the
        #                 same order as specified in the rule.  Contains
        #                 the filename that contains the data for this
        #                 object.  Also the attributes can be useful.
        # out_attributes  A dictionary with the desired attributes of the
        #                 output object.  Used so that the function knows
        #                 what to do, e.g. what algorithm to use for
        #                 normalization.
        # user_options    A dictionary with optional attributes specified
        #                 by the user.  For example, which GSEID to
        #                 download.
        # num_cores       An integer for the maximum number of cores to
        #                 use.
        # outfile         Full path of file to save the data.
        raise AssertionError, "No run function provided."


    def name_outfile(self, antecedents, user_options):
        # Come up with a meaningful name for the output file.  Returns
        # a local filename or a directory if this module requires one.  
        # Should not be a full directory path.
        #
        # No default implementation.  The outfile name should be
        # meaningful for the module.
        #
        # antecedents     As above.  Is useful when the name of the file
        #                 depends on the attributes of the input data.
        # user_options    As above.  Is useful when the name of the file
        #                 depends on the options specified by the user,
        #                 e.g. the GEO ID.
        raise AssertionError, "No name_outfile function provided."

    def set_out_attributes(self, antecedents, out_attributes):
        # Return the out_attributes that describes the output data.  This
        # is needed when the module looks at the data and sets some value,
        # e.g. is_logged.
        # 
        # The default implementation does not change the out_attributes.
        #
        # antecedents     As above.  Provides the filename for the input
        #                 data.
        # out_attributes  As above.
        return out_attributes


    ## def find_antecedents(
    ##     self, network, module_id, out_attributes, user_attributes, pool):
    ##     # Pull out the DataObjects from the pool that are the antecedents
    ##     # to module_id.  Returns a single or list of DataObjects.
    ##     #
    ##     # The default implementation can handle ModuleNodes with 1
    ##     # input data, or those where each input data has a different
    ##     # type.  If there are multiple inputs with different types,
    ##     # then you need to override this.
    ##     #
    ##     # network         As above.
    ##     # module_id       ID of module to get antecedents for.
    ##     # out_attributes  As above.  Useful when a module might have
    ##     #                 several antecedents of the same DataType, but
    ##     #                 different by attributes.  Can try to find one
    ##     #                 that fits the desired output attribute.  For
    ##     #                 example, different kinds of processing might
    ##     #                 lead to PCA plot.  But we'll still need to be
    ##     #                 able to distinguish them later.
    ##     # user_attributes As above.  Not ever used.
    ##     # pool            List of module_id -> DataObject

    ##     # TODO: Implement a default that works for modules with multiple
    ##     # antecedents.
    ##     from Betsy import module_utils

    ##     node = network.nodes[module_id]
    ##     # Make sure we don't have multiple input nodes with the same
    ##     # type.
    ##     seen = {}
    ##     for dt in node.in_datatypes:
    ##         assert dt.name not in seen, "Need to implement."
    ##         seen[dt.name] = 1

    ##     # Make sure each antecedent has the right data type.
    ##     filters = []
    ##     if len(node.in_datatypes) > 1:
    ##         for dt in node.in_datatypes:
    ##             x = module_utils.AntecedentFilter(datatype_name=dt.name)
    ##             filters.append(x)
    ##     assert len(node.in_datatypes) == 1 or \
    ##            len(node.in_datatypes) == len(filters)
        
    ##     x = module_utils.find_antecedents(
    ##         network, module_id, user_attributes, pool, *filters)
    ##     return x
