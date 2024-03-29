from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_filename):
        import itertools
        from genomicode import SimpleVariantMatrix
        from genomicode import AnnotationMatrix
        from Betsy import module_utils as mlib

        summary_file = in_data.identifier
        metadata = {}

        #x = mlib.get_user_option(
        #    user_options, "nonsynonymous_and_stopgain_only",
        #    allowed_values=["no", "yes"])
        #nonsynonymous_and_stopgain_only = (x == "yes")

        min_alt_reads = mlib.get_user_option(
            user_options, "filter_by_min_alt_reads", not_empty=True,
            type=int)
        assert min_alt_reads >= 0 and min_alt_reads < 10000

        min_total_reads = mlib.get_user_option(
            user_options, "filter_by_min_total_reads", not_empty=True,
            type=int)
        assert min_total_reads >= 0 and min_total_reads < 10000
        
        min_vaf = mlib.get_user_option(
            user_options, "filter_by_min_vaf", not_empty=True,
            type=float)
        assert min_vaf >= 0.0 and min_vaf < 1.0

        #min_gq = mlib.get_user_option(
        #    user_options, "filter_by_min_GQ", not_empty=True, type=float)
        #assert min_gq >= 0 and min_gq < 1000

        assert min_total_reads or min_alt_reads, "No filter"

        matrix = SimpleVariantMatrix.read_as_am(summary_file)
        #var_matrix = SimpleVariantMatrix.read(summary_file)
        #call_matrix = var_matrix.call_matrix
        #annot_matrix = var_matrix.annot_matrix

        #annovar_matrix = None
        #for (name, matrix) in var_matrix.named_matrices:
        #    if "ExonicFunc.refGene" in matrix.headers:
        #        annovar_matrix = matrix
        #        break
        #assert annovar_matrix, "Missing annotation: ExonicFunc.refGene"

        # copy.deepcopy is very slow.  Try to avoid it.
        # Strategy:
        # 1.  Make a list of the changes to be made.
        # 2.  Save the filtered rows.
        # 3.  Make the changes.
        # 4.  Save the non-filtered rows.
        I_remove = {}     # i -> 1
        call_remove = {}  # i -> (sample, caller) -> 1

        #CHROM = matrix.header2annots["______Chrom"]
        #POS = matrix.header2annots["______Pos"]
        #POS = [int(x) for x in POS]
        #REF = matrix.header2annots["______Ref"]
        #ALT = matrix.header2annots["______Alt"]

        # Optimization: normalize the headers for the samples and callers.
        sc2header = {}  # (sample, caller) -> header_h
        for sc in itertools.product(matrix.samples, matrix.callers):
            sample, caller = sc
            header = "%s___%s___Ref/Alt/VAF" % (sample, caller)
            header_h = matrix.normalize_header(header)
            assert header_h
            sc2header[sc] = header_h
        
        for i in range(matrix.num_annots()):
            has_calls = False   # whether this row has any calls.
            for sc in itertools.product(matrix.samples, matrix.callers):
                sample, caller = sc

                header_h = sc2header[sc]
                call_str = matrix.header2annots[header_h][i]
                if not call_str:
                    continue
                call = SimpleVariantMatrix._parse_call(call_str)

                filt = False
                # filter_by_min_alt_reads
                if min_alt_reads > 0 and \
                   (call.num_alt is None or call.num_alt < min_alt_reads):
                    filt = True
                # filter_by_min_total_reads
                if min_total_reads > 0 and (
                    call.total is None or call.total < min_total_reads):
                    filt = True

                # filter_by_min_vaf
                if min_vaf >= 1E-6 and (
                    call.vaf is None or call.vaf < min_vaf):
                    filt = True

                if filt:
                    if i not in call_remove:
                        call_remove[i] = {}
                    call_remove[i][sc] = 1
                else:
                    has_calls = True

            # If this coordinate has no more calls, then remove the
            # whole row.
            if not has_calls:
                I_remove[i] = 1
        I_remove = sorted(I_remove)

        # Write out a matrix of the discarded rows.
        filtered_matrix = AnnotationMatrix.rowslice(matrix, I_remove)
        SimpleVariantMatrix.write_from_am("discarded.txt", filtered_matrix)
        
        # Remove the calls.
        for i in call_remove:
            for sc in call_remove[i]:
                header_h = sc2header[sc]
                call_str = matrix.header2annots[header_h][i]
                assert call_str
                matrix.header2annots[header_h][i] = ""
                
        # Which rows to keep.
        I_remove_dict = {}.fromkeys(I_remove)
        I_keep = [
            i for i in range(matrix.num_annots()) if i not in I_remove_dict]
        filtered_matrix = AnnotationMatrix.rowslice(matrix, I_keep)
        SimpleVariantMatrix.write_from_am(out_filename, filtered_matrix)

        ## ## Filter out synonymous variants.
        ## #if nonsynonymous_and_stopgain_only:
        ## #    # Make sure annotated with Annovar.
        ## #    assert "ExonicFunc.refGene" in annovar_matrix.headers
        ## #    exonic_func = annovar_matrix["ExonicFunc.refGene"]
        ## #    for i, efunc in enumerate(exonic_func):
        ## #        efunc = exonic_func[i]
        ## #        assert efunc in [
        ## #            "", "nonsynonymous SNV", "synonymous SNV",
        ## #            "stopgain", "stoploss",
        ## #            "frameshift substitution", "nonframeshift substitution",
        ## #            "unknown"], \
        ## #            "Unknown exonic_func: %s" % efunc
        ## #        if efunc not in ["nonsynonymous SNV", "stopgain"]:
        ## #            I_remove[i] = 1
        ## #            continue
                
        ## # Filter based on the calls.
        ## if min_alt_reads > 0 or min_total_reads > 0:
        ##     all_coord = call_matrix.coord2samplecaller2call.keys()
        ##     for coord in all_coord:
        ##         all_sc = call_matrix.coord2samplecaller2call[coord].keys()
        ##         for sc in all_sc:
        ##             # SimpleVariantMatrix.Call object.
        ##             call = call_matrix.coord2samplecaller2call[coord][sc]
                    
        ##             # filter_by_min_alt_reads
        ##             if min_alt_reads > 0 and \
        ##                (call.num_alt is None or call.num_alt < min_alt_reads):
        ##                 if coord not in call_remove:
        ##                     call_remove[coord] = {}
        ##                 call_remove[coord][sc] = 1
                        
        ##             # filter_by_min_total_reads
        ##             if min_total_reads > 0 and (
        ##                 call.total is None or call.total < min_total_reads):
        ##                 if coord not in call_remove:
        ##                     call_remove[coord] = {}
        ##                 call_remove[coord][sc] = 1
                        
        ## # Filter based on VAF.
        ## if min_vaf >= 1E-6:
        ##     all_coord = call_matrix.coord2samplecaller2call.keys()
        ##     for coord in all_coord:
        ##         all_sc = call_matrix.coord2samplecaller2call[coord].keys()
        ##         for sc in all_sc:
        ##             call = call_matrix.coord2samplecaller2call[coord][sc]

        ##             # filter_by_min_vaf
        ##             if call.vaf is None or call.vaf < min_vaf:
        ##                 if coord not in call_remove:
        ##                     call_remove[coord] = {}
        ##                 call_remove[coord][sc] = 1

        ## # If any of these coordinates have no more variants, then
        ## # remove the whole row.
        ## if call_remove:
        ##     chrom, pos = annot_matrix["Chrom"], annot_matrix["Pos"]
        ##     ref, alt = annot_matrix["Ref"], annot_matrix["Alt"]
        ##     pos = [int(x) for x in pos]
        ##     coord2i = {}
        ##     for i, coord in enumerate(zip(chrom, pos, ref, alt)):
        ##         coord2i[coord] = i

        ##     for coord in call_remove:
        ##         num_remove = len(call_remove[coord])
        ##         num_calls = len(call_matrix.coord2samplecaller2call[coord])
        ##         assert num_remove <= num_calls
        ##         if num_remove == num_calls:
        ##             i = coord2i[coord]
        ##             I_remove[i] = 1

        ## # Make a matrix of the discarded rows.
        ## old_annot_matrix = var_matrix.annot_matrix
        ## old_named_matrices = var_matrix.named_matrices
        ## filtered_matrix = var_matrix
        ## x = AnnotationMatrix.rowslice(var_matrix.annot_matrix, I_remove)
        ## filtered_matrix.annot_matrix = x
        ## named_matrices = []
        ## for (name, matrix) in var_matrix.named_matrices:
        ##     matrix = AnnotationMatrix.rowslice(matrix, I_remove)
        ##     named_matrices.append((name, matrix))
        ## filtered_matrix.named_matrices = named_matrices
        ## SimpleVariantMatrix.write("discarded.txt", filtered_matrix)
        ## var_matrix.annot_matrix = old_annot_matrix
        ## var_matrix.named_matrices = old_named_matrices
        
        ## # Remove the calls.
        ## for coord in call_remove:
        ##     chrom, pos, ref, alt = coord
        ##     for (sample, caller) in call_remove[coord]:
        ##         var_matrix.call_matrix.set_call(
        ##             chrom, pos, ref, alt, sample, caller, None)

        ## # Which rows to keep.
        ## I_keep = [
        ##     i for i in range(var_matrix.num_variants()) if i not in I_remove]
        ## # Filter annotation matrix
        ## var_matrix.annot_matrix = AnnotationMatrix.rowslice(
        ##     var_matrix.annot_matrix, I_keep)
        ## # Filter named matrices.
        ## for i, (name, matrix) in enumerate(var_matrix.named_matrices):
        ##     matrix = AnnotationMatrix.rowslice(matrix, I_keep)
        ##     var_matrix.named_matrices[i] = (name, matrix)
        
        ## SimpleVariantMatrix.write(out_filename, var_matrix)

        return metadata
            
    
    def name_outfile(self, antecedents, user_options):
        return "calls.txt"
