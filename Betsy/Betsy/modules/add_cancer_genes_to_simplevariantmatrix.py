from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_filename):
        from genomicode import filelib
        from genomicode import AnnotationMatrix
        from genomicode import SimpleVariantMatrix
        from Betsy import module_utils as mlib

        simple_node = in_data
        filelib.assert_exists_nz(simple_node.identifier)

        gene_file = mlib.get_user_option(
            user_options, "cancer_genes_file", not_empty=True, check_file=True)

        # Read the cancer genes file.
        # <Gene ID>  <Gene Symbol>  <Dataset>  ...
        symbol2info = {}  # symbol -> d
        gene_iter = filelib.read_row(gene_file, header=1)
        for d in gene_iter:
            assert "Gene Symbol" in d._header
            if not d.Gene_Symbol:
                continue
            symbol2info[d.Gene_Symbol] = d

        # Read the variant file.
        SVM = SimpleVariantMatrix.read(simple_node.identifier)

        GENE_H = "Gene.refGene"
        annovar_matrix = None
        for (name, matrix) in SVM.named_matrices:
            if GENE_H in matrix.headers:
                annovar_matrix = matrix
                break
        assert annovar_matrix, "Missing annotation: %s" % GENE_H
        GENES = annovar_matrix[GENE_H]

        # Align the matrix to the simple variant matrix.
        headers = gene_iter._header[2:]
        assert headers
        matrix = [[None]*len(headers) for i in range(len(GENES))]
        for i, gene_str in enumerate(GENES):
            # Format of genes:     Format of output
            # PFN1P2                  5.2
            # PMS2P2,PMS2P7           2.2,8.6
            # If a gene is missing, then skip it.

            values = [""] * len(headers)
            # If any of these genes are in a gene set, then put a "1"
            # in the matrix.
            genes = gene_str.split(",")
            for gene in genes:
                if gene not in symbol2info:
                    continue
                d = symbol2info[gene]
                info = d._cols[2:]
                assert len(info) == len(values)
                for j in range(len(info)):
                    assert info[j] in ["", "1"]
                    if info[j] == "1":
                        values[j] = 1
            matrix[i] = values

        # Add the matrix back to the simple variant matrix.
        all_annots = []
        for j in range(len(headers)):
            x = [matrix[i][j] for i in range(len(matrix))]
            all_annots.append(x)
        x = AnnotationMatrix.create_from_annotations(headers, all_annots)
        SVM.named_matrices.append(("Cancer Genes", x))

        # Write to file.
        SimpleVariantMatrix.write(out_filename, SVM)
                    
                    
    def name_outfile(self, antecedents, user_options):
        return "matrix.txt"
