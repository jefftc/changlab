from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_filename):
        import arrayio
        from genomicode import filelib
        from genomicode import AnnotationMatrix
        from genomicode import SimpleVariantMatrix

        simple_node, signal_node = antecedents
        filelib.assert_exists_nz(simple_node.identifier)
        filelib.assert_exists_nz(signal_node.identifier)

        # Read the variant file.
        SVM = SimpleVariantMatrix.read(simple_node.identifier)
        #AM = SVM.annot_matrix
        #assert GENE_H in AM.headers

        # Read the gene expression file.
        GXP = arrayio.read(signal_node.identifier)

        # Make sure the samples from the variant matrix can be found
        # in the gene expression matrix.
        GXP_samples = GXP.col_names(arrayio.COL_ID)
        missing = [x for x in SVM.samples if x not in GXP_samples]
        assert len(missing) < len(SVM.samples), (
            "SimpleVariantMatrix and gene expression file have "
            "no common samples.")
        # Actually, may not have all the same samples.  For example, a
        # gene expression profile might not have been calculated for
        # the germline sample.  So ignore if something is missing.
        #x = missing
        #if len(x) > 5:
        #    x = x[:5] + ["..."]
        #msg = "Samples (%d) not found in gene expression file: %s" % (
        #    len(missing), ", ".join(x))
        #assert not missing, msg

        # Add all the samples from the gene expression file.
        SAMPLES = GXP_samples

        # Find the genes in each row.
        GENE_H = "Gene.refGene"
        annovar_matrix = None
        for (name, matrix) in SVM.named_matrices:
            if GENE_H in matrix.headers:
                annovar_matrix = matrix
                break
        assert annovar_matrix, "Missing annotation: %s" % GENE_H
        GENES = annovar_matrix[GENE_H]

        # Make a list of the genes.
        genes = {}
        for i, gene_str in enumerate(GENES):
            # Format of genes:
            # PFN1P2
            # PMS2P2,PMS2P7
            for x in gene_str.split(","):
                genes[x] = 1
        genes = sorted(genes)

        # Make a matrix of the gene expression values for each gene
        # and each sample.
        #I = [GXP_samples.index(x) for x in SVM.samples]
        #GXP_a = GXP.matrix(genes, I)  # align the matrices.
        GXP_a = GXP.matrix(genes, None)
        
        # Write out the expression matrix for debugging purposes.
        arrayio.write(GXP_a, "expression.txt")

        # Search for each of the genes in the matrix.
        gene2I = {}   # gene -> list of row indexes
        for gene in genes:
            x = GXP_a._index(row=gene)
            I_row, i_col = x
            if I_row:
                gene2I[gene] = I_row

        # Align the gene expression matrix to the simple variant
        # matrix.
        #matrix = [[None]*len(SVM.samples) for i in range(len(GENES))]
        matrix = [[None]*len(SAMPLES) for i in range(len(GENES))]
        for i, gene_str in enumerate(GENES):
            # Format of genes:     Format of output
            # PFN1P2                  5.2
            # PMS2P2,PMS2P7           2.2,8.6
            # If a gene is missing, then skip it.
            genes = gene_str.split(",")
            #for j in range(len(SVM.samples)):
            for j in range(len(SAMPLES)):
                values = []  # expression values for each gene.
                for k in range(len(genes)):
                    if genes[k] not in gene2I:
                        continue
                    x = [GXP_a._X[l][j] for l in gene2I[genes[k]]]
                    # If there are multiple instances of this gene,
                    # then pick the one with the maximum expression.
                    x = max(x)
                    values.append(x)
                values = [_pretty_gxp(x) for x in values]
                x = ",".join(values)
                matrix[i][j] = x

        # Add the matrix back to the simple variant matrix.
        #headers = SVM.samples
        headers = SAMPLES
        all_annots = []
        for j in range(len(headers)):
            x = [matrix[i][j] for i in range(len(matrix))]
            all_annots.append(x)
        x = AnnotationMatrix.create_from_annotations(headers, all_annots)
        SVM.named_matrices.append(("Gene Expression", x))

        # Write to file.
        SimpleVariantMatrix.write(out_filename, SVM)

        
        ## x = ["%s Exp" % x for x in SVM.samples]
        ## headers = AM.headers + x
        ## all_annots = []
        ## for h in AM.headers_h:
        ##     all_annots.append(AM.header2annots[h])
        ## for i in range(len(SVM.samples)):
        ##     x = [x[i] for x in matrix]
        ##     assert len(x) == AM.num_annots()
        ##     all_annots.append(x)
        ## assert not AM.headerlines
        ## x = AnnotationMatrix.create_from_annotations(headers, all_annots)
        ## x = SimpleVariantMatrix.SimpleVariantMatrix(
        ##     SVM.samples, SVM.callers, x, SVM.call_matrix)
        ## SimpleVariantMatrix.write(out_filename, x)
                    
                    
    def name_outfile(self, antecedents, user_options):
        return "matrix.txt"


def _pretty_gxp(value):
    # Format a gene expression value.
    if value < 1E-5:  # Don't have values like 0.0.
        return "0"
    return "%.3f" % value
    
