from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        import arrayio
        #from genomicode import config
        #from genomicode import arrayannot
        from genomicode import arrayplatformlib
        from genomicode import parallel
        from genomicode import filelib
        from Betsy import module_utils as mlib

        M = arrayio.read(in_data.identifier)
        metadata = {}

        # Add GENE_ID, GENE_SYMBOL, and DESCRIPTION.  Figure out which
        # platforms provide each one of this.
        CATEGORIES = [
            arrayplatformlib.GENE_ID,
            arrayplatformlib.GENE_SYMBOL,
            # biomaRt doesn't convert description.  So just ignore it
            # for now.
            # TODO: implement DESCRIPTION.
            #arrayplatformlib.DESCRIPTION,
            ]

        #all_platforms = arrayplatformlib.identify_all_platforms_of_matrix(M)
        #assert all_platforms, "Unknown platform: %s" % in_data.identifier
        #header, platform_name = all_platforms[0]
        scores = arrayplatformlib.score_matrix(M)
        assert scores, "I could not identify platform."
        score = scores[0]
        assert score.max_score >= 0.75, "I could not identify platform."
        platform = arrayplatformlib.find_platform_by_name(score.platform_name)
        
        to_add = []  # list of platform names
        for category in CATEGORIES:
            x = arrayplatformlib.PLATFORMS
            x = [x for x in x if x.category == category]
            x = [x for x in x if x.bm_organism == platform.bm_organism]
            x = [x for x in x if x.name != score.platform_name]
            # Take the first one, of any.
            if x:
                to_add.append(x[0].name)
        assert to_add, "No platforms to add"

        annotate = mlib.get_config("annotate_matrix", which_assert_file=True)
        sq = parallel.quote
        cmd = [
            "python",
            sq(annotate),
            "--header", sq(score.header),
            ]
        for x in to_add:
            x = ["--platform", sq(x)]
            cmd.extend(x)
        cmd.append(in_data.identifier)
        cmd = " ".join(cmd)
        x = parallel.sshell(cmd)
        metadata["commands"] = [cmd]
        assert x.find("Traceback") < 0, x
        
        open(outfile, 'w').write(x)

        ## ids = M._row_order
        ## probe_header = all_platforms[0][0]
        ## probe_id = M._row_names[probe_header]
        ## new_ids = ids[:]
        ## new_ids.remove(probe_header)
        ## #annotate_type = parameters['annotate_type']
        ## #if annotate_type == 'all':
        ## annotate_header = arrayplatformlib.annotate_header
        ## #elif annotate_type == 'gene_id':
        ## #    annotate_header = ['Gene ID']
        ## dictionary = arrayannot.annotate_probes_multiple(
        ##     probe_id, annotate_header)
        
        ## column = []
        ## for id_ in new_ids:
        ##     flag = True
        ##     for key in dictionary.keys():
        ##         flag = True
        ##         value_list = dictionary[key]
        ##         gene_list = M._row_names[id_]
        ##         for gene in gene_list:
        ##             if gene not in value_list:
        ##                 flag = False
        ##                 break
        ##         if flag:
        ##             column.append((key, id_))
    
        ## header = [i[0] for i in column]
        ## miss_header = list(set(annotate_header).difference(set(header)))
        ## original_ids = ids[:]
        ## for col in miss_header:
        ##     col_2 = col
        ##     if col in original_ids:
        ##         col_1 = col + '_1'
        ##         col_2 = col + '_2'
        ##         M = module_utils.replace_matrix_header(M, col, col_1)
        ##         ids = M._row_order
        ##     ids.append(col_2)
        ##     M._row_order = ids
        ##     M._row_names[col_2] = dictionary[col]
        
        ## f = file(outfile, 'w')
        ## arrayio.tab_delimited_format.write(M, f)
        ## f.close()

        return metadata
        

    def name_outfile(self, antecedents, user_options):
        return "signal_annot.tdf"
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #filename = 'signal_annot_' + original_file + '.tdf'
        #return filename
