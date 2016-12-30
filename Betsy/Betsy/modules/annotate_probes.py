from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        import StringIO
        import arrayio
        from genomicode import arrayplatformlib
        from genomicode import parallel
        from genomicode import filelib
        from genomicode import AnnotationMatrix
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
            "--no_na", 
            "--header", sq(score.header),
            ]
        for x in to_add:
            x = ["--platform", sq(x)]
            cmd.extend(x)
        cmd.append(in_data.identifier)
        cmd = " ".join(cmd)
        data = parallel.sshell(cmd)
        metadata["commands"] = [cmd]
        assert data.find("Traceback") < 0, data

        # Clean up the headers.
        platform2pretty = {
            "Entrez_ID_human" : "Gene ID",
            "Entrez_Symbol_human" : "Gene Symbol",
            "Entrez_ID_mouse" : "Gene ID",
            "Entrez_Symbol_mouse" : "Gene Symbol",
            }
        handle = open(outfile, 'w')
        header_written = False
        for cols in filelib.read_cols(StringIO.StringIO(data)):
            if not header_written:
                cols = [platform2pretty.get(x, x) for x in cols]
                cols = AnnotationMatrix.uniquify_headers(cols)
                header_written = True
            print >>handle, "\t".join(cols)

        return metadata
        

    def name_outfile(self, antecedents, user_options):
        return "signal_annot.tdf"
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #filename = 'signal_annot_' + original_file + '.tdf'
        #return filename
