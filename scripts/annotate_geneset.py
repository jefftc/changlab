#!/usr/bin/env python

# Functions:
# read_geneset
# read_annotations


import os


def read_geneset(geneset):
    # Parse a geneset specified by the user.  geneset is in the format
    # of <filename>[,<geneset>,<geneset>,...].  Return a list of
    # (<filename>, <geneset>, list of genes).
    from genomicode import genesetlib
    
    x = geneset.split(",")
    assert len(x) >= 1
    filename, genesets = x[0], x[1:]
    assert os.path.exists(filename), "I could not find file: %s" % filename
    
    data = []
    for geneset in genesets:
        genes = genesetlib.read_genes(filename, geneset, allow_tdf=True)
        assert genes, "I could not find the gene set: %s" % geneset
        x = filename, geneset, genes
        data.append(x)
    return data


def read_annotations(annotation_file):
    # Parse the annotations specified by the user.  annotation_file is
    # a file of gene sets.  Return a dictionary where the keys are the
    # names of the annotations and the values are a list of genes.
    from genomicode import genesetlib

    assert os.path.exists(annotation_file), \
           "I could not find file: %s" % annotation_file

    annot2genes = {}
    for x in genesetlib.read_genesets(annotation_file):
        name, description, genes = x
        annot2genes[name] = genes
    assert annot2genes, "No annotations."
    return annot2genes


def read_gene_descriptor(gene_descriptor, geneset):
    # Read pretty names for the genes.  gene_descriptor is in the
    # format of <filename>,<header>.  Return a dictionary of gene ->
    # pretty name.
    
    from genomicode import genesetlib

    x = gene_descriptor.split(",")
    assert len(x) >= 2
    filename, pretty_header = x
    assert os.path.exists(filename), "I could not find file: %s" % filename

    header2genes = {}
    for x in genesetlib.read_tdf(
        filename, preserve_spaces=True, allow_duplicates=True):
        name, description, genes = x
        header2genes[name] = genes
    assert header2genes, "No genes."

    # Find the column that contains the genes in the gene set.  Since
    # some of the genes may not be annotated, provide some leeway
    # here.
    header2counts = {}
    for (header, genes) in header2genes.iteritems():
        count = len(set(geneset).intersection(genes))
        header2counts[header] = count
    best_header = best_count = None
    for (header, count) in header2counts.iteritems():
        if best_count is None or count > best_count:
            best_header, best_count = header, count
    assert best_count >= len(geneset)/2.0, \
                  "I could not find the genes in the descriptor file."

    gene_header = best_header
    genes = header2genes[gene_header]
    pretty = header2genes[pretty_header]
    assert len(genes) == len(pretty)

    gene2pretty = {}
    for g, p in zip(genes, pretty):
        if not g or not p:
            continue
        gene2pretty[g] = p
    
    return gene2pretty


def read_annotation_descriptor(annotation_descriptor, annotations):
    # Return list of (header, annots), where the annots are aligned to
    # the annotations.
    from genomicode import jmath
    from genomicode import genesetlib
    
    filename = annotation_descriptor
    assert os.path.exists(filename), "I could not find file: %s" % filename

    header2annots = []  # list of (header, annots)
    for x in genesetlib.read_tdf(
        filename, preserve_spaces=True, allow_duplicates=True):
        name, description, annots = x
        header2annots.append((name, annots))
    assert header2annots, "No annots."

    # Find the column that contains the annotations.
    header2counts = {}
    for (header, annots) in header2annots:
        count = len(set(annotations).intersection(annots))
        header2counts[header] = count
    best_header = best_count = None
    for (header, count) in header2counts.iteritems():
        if best_count is None or count > best_count:
            best_header, best_count = header, count
    assert best_count >= len(annotations)/2.0, \
                  "I could not find the annotations in the descriptor file."
    annot_header = best_header

    # Align the annotation matrix to the annotations.
    annot_annots = None
    for (header, annots) in header2annots:
        if header == annot_header:
            annot_annots = annots
    assert annot_annots
    I = jmath.match(annotations, annot_annots)
    
    header2annots_aligned = []
    for header, annots in header2annots:
        annots_aligned = []
        for i in I:
            if i is None:
                annots_aligned.append("")
            else:
                annots_aligned.append(annots[i])
        x = header, annots_aligned
        header2annots_aligned.append(x)
    header2annots = header2annots_aligned

    return header2annots


def _annotate_genes_h(
    geneset, background, annot2genes, annotations, min_genes):
    import math
    from genomicode import jmath

    scores = {}
    for annot in annotations:
        # Make a list of all the genes with this annotation (that
        # shows up in our background list.
        x = annot2genes[annot]
        x = [x for x in x if x in background]
        annot_genes = x
        if not annot_genes:
            continue
        if min_genes is not None and len(annot_genes) < min_genes:
            continue

        total_genes = len(background)
        genes_in_list = len(geneset)
        genes_with_annot = len(annot_genes)

        # L 0/1  Gene in list
        # A 0/1  Gene with annotation
        L1A1 = len([x for x in geneset if x in annot_genes])
        L1A0 = genes_in_list - L1A1
        L0A1 = genes_with_annot - L1A1
        L0A0 = total_genes - L1A1 - L1A0 - L0A1
        assert L1A1 >= 0
        assert L1A0 >= 0
        assert L0A1 >= 0
        assert L0A0 >= 0
        perc_list = float(L1A1) / (L1A1+L1A0)
        perc_background = float(L0A1) / (L0A1+L0A0)
        if not perc_background:
            perc_list = (float(L1A1)+1) / (L1A1+L1A0+1)
            perc_background = (float(L0A1)+1) / (L0A1+L0A0+1)
        fe = perc_list / perc_background

        p_value = jmath.fisher_test(L0A0, L0A1, L1A0, L1A1)
        nl10p = -math.log(p_value, 10)

        x = [x for x in geneset if x in annot_genes]
        x = sorted(x)
        L1A1_genes = x

        
        x = L1A1, L1A0, L0A1, L0A0, fe, nl10p, L1A1_genes
        scores[annot] = x
    return scores


def annotate_genes(geneset, background, annot2genes, min_genes, num_procs):
    import multiprocessing
    
    pool = multiprocessing.Pool(num_procs)

    all_annotations = sorted(annot2genes)
    proc2annots = {}
    for i, annot in enumerate(all_annotations):
        proc = i % num_procs
        if proc not in proc2annots:
            proc2annots[proc] = []
        proc2annots[proc].append(annot)

    scores = {}
    results = []
    for annots in proc2annots.itervalues():
        fn = _annotate_genes_h
        args = geneset, background, annot2genes, annots, min_genes
        keywds = {}
        if num_procs == 1:
            x = fn(*args)
            scores.update(x)
        else:
            x = pool.apply_async(fn, args, keywds)
            results.append(x)
    pool.close()
    pool.join()
    for x in results:
        x = x.get()
        scores.update(x)

    return scores

    
def main():
    import argparse
    from genomicode import jmath

    parser = argparse.ArgumentParser(
        description="Annotate a gene set with Gene Ontology codes.")
    parser.add_argument(
        "-j", dest="num_procs", type=int, default=1,
        help="Number of jobs to run in parallel.")
    parser.add_argument(
        "--min_genes", type=int, default=0,
        help="Ignore annotations that do not have at least this number of "
        "genes.")

    group = parser.add_argument_group(title="Required arguments")
    group.add_argument(
        "--geneset", default=None,
        help="Annotate this geneset.  If multiple gene sets are provided, "
        "their genes will be combined.  "
        "Format: <gmx/gmt_file>[,<geneset>,<geneset>,...]")
    group.add_argument(
        "--background_geneset", default=None,
        help="Genes are selected from this background geneset.  "
        "Format: <gmx/gmt_file>,<geneset>")
    group.add_argument(
        "--annotation", default=None,
        help="A gene set file that contain the annotations for the genes.  "
        "Format: <gmx/gmt_file>,<geneset>")

    group = parser.add_argument_group(title="Descriptors")
    group.add_argument(
        "--annotation_descriptor", default=None,
        help="A text file that contains more information about the "
        "annotations.  "
        "One of the columns should contain descriptors that match the "
        "names of the annotations given in the annotation file.  ")
    group.add_argument(
        "--gene_descriptor", default=None,
        help="A text file that contains alternate names for the genes.  "
        "One of the columns should contain IDs that match geneset.  "
        "gene_name_header should contain the name of the genes to show "
        "in the output file.  "
        "Format: <filename>,<gene_name_header>")
    
    
    args = parser.parse_args()
    assert args.geneset, "Please specify a gene set."
    assert args.background_geneset, "Please specify a background gene set."
    assert args.annotation, "Please specify an annotation file."
    if args.num_procs < 1 or args.num_procs > 100:
        parser.error("Please specify between 1 and 100 processes.")

    # Read the gene sets.
    x1 = read_geneset(args.geneset)
    x2 = read_geneset(args.background_geneset)
    assert len(x1) >= 1
    assert len(x2) == 1
    background = x2[0][-1]
    # Save the gene sets individually.
    geneset2genes = {}
    for x in x1:
        filename, gs_name, genes = x
        geneset2genes[gs_name] = genes
    # Combine the genes in the gene sets.
    geneset = []
    for genes in geneset2genes.itervalues():
        geneset.extend(genes)
    geneset = sorted({}.fromkeys(geneset))

    # Make sure each gene is in the background.
    missing = [x for x in geneset if x not in background]
    if len(missing) == len(geneset):
        assert False, "All genes from the geneset are missing in background."
    elif missing and len(missing) <= 5:
        x = ", ".join(missing)
        assert False, "Genes are missing from the background: %s" % x
    elif missing:
        assert False, "%d genes are missing from the background." % len(
            missing)

    # Read the annotations.
    x = read_annotations(args.annotation)
    annot2genes = x

    # Calculate the score for each of the annotations.
    scores = annotate_genes(
        geneset, background, annot2genes, args.min_genes, args.num_procs)
    all_annots = sorted(scores)

    # Do multiple hypothesis correction.
    nl10ps = []
    for x in scores.itervalues():
        L1A1, L1A0, L0A1, L0A0, fe, nl10p, L1A1_genes = x
        nl10ps.append(nl10p)
    p_values = [10**-x for x in nl10ps]
    bonferroni = jmath.cmh_bonferroni(p_values)
    fdr = jmath.cmh_fdr_bh(p_values)
    nl10p2stats = {}  # nl10p -> p_value, bonf, fdr
    for n, p, b, f in zip(nl10ps, p_values, bonferroni, fdr):
        nl10p2stats[n] = p, b, f

    # Read the descriptors for output.
    gene2pretty = {}
    if args.gene_descriptor:
        gene2pretty = read_gene_descriptor(args.gene_descriptor, background)
    if args.annotation_descriptor:
        annot_descriptors = read_annotation_descriptor(
            args.annotation_descriptor, all_annots)

    annot_headers = [x[0] for x in annot_descriptors]
    all_genesets = sorted(geneset2genes)
    header = annot_headers + [
        "Annotation",
        "Your Genes (With Ann)", "Your Genes (No Ann)",
        "Other Genes (With Ann)", "Other Genes (No Ann)",
        "Your Genes Annotated", "Other Genes Annotated",
        "Fold Enrichment", "neg log_10(p value)", "Bonferroni", "FDR"] + \
        all_genesets
    print "\t".join(header)
    for i, annot in enumerate(all_annots):
        L1A1, L1A0, L0A1, L0A0, fe, nl10p, L1A1_genes = scores[annot]
        perc_geneset = float(L1A1) / (L1A1+L1A0)
        perc_background = float(L0A1) / (L0A1+L0A0)
        pvalue, bonf, fdr = nl10p2stats[nl10p]

        # Pull out the descriptors for the annotations.
        annot_info = [x[1][i] for x in annot_descriptors]

        # Separate the L1A1_genes into individual genes.
        geneset2L1A1 = {}
        for gs_name, genes in geneset2genes.iteritems():
            x = [x for x in L1A1_genes if x in genes]
            geneset2L1A1[gs_name] = x

        # Convert the gene names to pretty names.
        for gs_name, genes in geneset2L1A1.iteritems():
            x = [gene2pretty.get(x, x) for x in genes]
            if gene2pretty:
                x.sort()
            geneset2L1A1[gs_name] = x

        # Format the genes for output in the table.
        geneset_list = []
        for gs_name in all_genesets:
            x = geneset2L1A1.get(gs_name, [])
            x = " ".join(x)
            geneset_list.append(x)

        x = annot_info + [
            annot, L1A1, L1A0, L0A1, L0A0, perc_geneset, perc_background,
            fe, nl10p, bonf, fdr] + geneset_list
        assert len(x) == len(header), "%d %d" % (len(x), len(header))
        print "\t".join(map(str, x))

    

if __name__ == '__main__':
    #import cProfile as profile
    #profile.runctx("main()", globals(), locals())
    main()
