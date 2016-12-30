from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import filelib
        from genomicode import hashlib
        from genomicode import jmath
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
        header = None
        for d in gene_iter:
            assert "Gene Symbol" in d._header
            if header is None:
                header = [
                    x for x in d._header
                    if x not in ["Gene ID", "Gene Symbol"]]
            if not d.Gene_Symbol:
                continue
            symbol2info[d.Gene_Symbol] = d

        # Read the variant file.
        SVM = SimpleVariantMatrix.read_as_am(simple_node.identifier)

        GENE_H = "Annovar______Gene.refGene"
        assert GENE_H in SVM.headers, "Missing annotation: %s" % GENE_H
        GENES = SVM[GENE_H]

        # Align the matrix to the simple variant matrix.
        gene_headers = header
        gene_annotations = []
        for i, gene_str in enumerate(GENES):
            # Format of genes:
            # PFN1P2
            # PMS2P2,PMS2P7
            values = [""] * len(gene_headers)
            genes = gene_str.split(",")
            for gene in genes:
                if gene not in symbol2info:
                    continue
                d = symbol2info[gene]
                for j, h in enumerate(gene_headers):
                    h = hashlib.hash_var(h)
                    assert hasattr(d, h)
                    x = getattr(d, h)
                    assert x in ["", "1"]
                    if x == "1":
                        values[j] = 1
            gene_annotations.append(values)
        # Convert the headers and annotations to SVM format.
        gene_headers = ["Cancer Genes______%s" % x for x in gene_headers]
        gene_annotations = jmath.transpose(gene_annotations)

        # Make the new SimpleVariantMatrix.
        # Figure out where to put these annotations.
        INDEX = 4
        # If Annovar exists, put after.
        I = [i for (i, x) in enumerate(SVM.headers)
             if x.upper().startswith("ANNOVAR")]
        if I:
            INDEX = max(INDEX, max(I)+1)
        # If SnpEff exists, put after.
        I = [i for (i, x) in enumerate(SVM.headers)
             if x.upper().startswith("SNPEFF")]
        if I:
            INDEX = max(INDEX, max(I)+1)
        # If COSMIC exists, put after.
        I = [i for (i, x) in enumerate(SVM.headers)
             if x.upper().startswith("COSMIC")]
        if I:
            INDEX = max(INDEX, max(I)+1)
        headers = SVM.headers[:INDEX] + gene_headers + SVM.headers[INDEX:]
        x = [SVM.header2annots[x] for x in SVM.headers_h]
        all_annots = x[:INDEX] + gene_annotations + x[INDEX:]
        merged = AnnotationMatrix.create_from_annotations(
            headers, all_annots, headerlines=SVM.headerlines)

        SimpleVariantMatrix.write_from_am(outfile, merged)
        

        
##         headers = gene_iter._header[2:]
##         assert headers
##         matrix = [[None]*len(headers) for i in range(len(GENES))]
##         for i, gene_str in enumerate(GENES):
##             # Format of genes:     Format of output
##             # PFN1P2                  5.2
##             # PMS2P2,PMS2P7           2.2,8.6
##             # If a gene is missing, then skip it.

##             values = [""] * len(headers)
##             # If any of these genes are in a gene set, then put a "1"
##             # in the matrix.
##             genes = gene_str.split(",")
##             for gene in genes:
##                 if gene not in symbol2info:
##                     continue
##                 d = symbol2info[gene]
##                 info = d._cols[2:]
##                 assert len(info) == len(values)
##                 for j in range(len(info)):
##                     assert info[j] in ["", "1"]
##                     if info[j] == "1":
##                         values[j] = 1
##             matrix[i] = values

##         # Add the matrix back to the simple variant matrix.
##         all_annots = []
##         for j in range(len(headers)):
##             x = [matrix[i][j] for i in range(len(matrix))]
##             all_annots.append(x)
##         x = AnnotationMatrix.create_from_annotations(headers, all_annots)
##         SVM.named_matrices.append(("Cancer Genes", x))

##         # Write to file.
##         SimpleVariantMatrix.write(out_filename, SVM)
                    
                    
    def name_outfile(self, antecedents, user_options):
        return "matrix.txt"
