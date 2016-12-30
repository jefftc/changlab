from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import os
        from genomicode import filelib
        from genomicode import parallel
        from Betsy import module_utils as mlib

        svm_node, vcf_node = antecedents
        vcf_filenames = filelib.list_files_in_path(
            vcf_node.identifier, endswith=".vcf", not_empty=True)
        metadata = {}

        # 1.  vcf_filenames
        # 2.  parsed_snpeff_files   one for each VCF file
        # 3.  merged_snpeff_file    just one file
        # 4.  clean_snpeff_file     clean up the annotations to final form
        # 5.  outfile

        merged_snpeff_file = "snpeff.merged.txt"
        cleaned_snpeff_file = "snpeff.clean.txt"

        jobs = []
        for vcf_filename in vcf_filenames:
            path, caller, ext = mlib.splitpath(vcf_filename)
            parsed_snpeff_file = "%s.parsed.txt" % caller
            j = filelib.GenericObject(
                caller=caller,
                vcf_filename=vcf_filename,
                parsed_snpeff_file=parsed_snpeff_file,
                )
            jobs.append(j)

        # Parse each of the snpeff files.
        commands = []
        for j in jobs:
            args = j.vcf_filename, j.parsed_snpeff_file
            # Debugging.  If this file exists, do not generate it
            # again.
            if os.path.exists(j.parsed_snpeff_file):
                continue
            x = parse_snpeff_file, args, {}
            commands.append(x)
        parallel.pyfun(commands, num_procs=num_cores)
        metadata["num_cores"] = num_cores

        # Merge the parsed files.
        x = [j.parsed_snpeff_file for j in jobs]
        x = [x for x in x if os.path.exists(x)]
        parsed_files = x
        # For debugging, don't regenerate if I don't need to.
        if not filelib.exists_nz(merged_snpeff_file):
            merge_parsed_files(parsed_files, merged_snpeff_file)

        # Clean up the snpEff file.  Coordinates should be unique.
        # For debugging, don't regenerate if I don't need to.
        if not filelib.exists_nz(cleaned_snpeff_file):
            clean_snpeff_file(merged_snpeff_file, cleaned_snpeff_file)

        # Merge the snpEff annotations into the SimpleVariantMatrix.
        add_snpeff_to_svm(svm_node.identifier, cleaned_snpeff_file, outfile)

        return metadata


    def name_outfile(self, antecedents, user_options):
        return "svm.snpeff_annotated.txt"


def parse_snpeff_file(vcf_filename, out_filename):
    from genomicode import vcflib

    # Parse out the snpEff annotations.  Should have ANN in INFO.
    # Make a tab-delimited text file containin columns:
    # Chrom  Pos  Ref  Alt  <snpEff-specific columns>
    #
    # ##INFO=<ID=ANN,Number=.,Type=String,
    #     Description="Functional annotations: '
    #     Allele |
    #     Annotation |
    #     Annotation_Impact |
    #     Gene_Name |
    #     Gene_ID |
    #     Feature_Type |
    #     Feature_ID |
    #     Transcript_BioType |
    #     Rank |
    #     HGVS.c |
    #     HGVS.p |
    #     cDNA.pos / cDNA.length |
    #     CDS.pos / CDS.length |
    #     AA.pos / AA.length |
    #     Distance | ERRORS / WARNINGS / INFO' ">

    vcf = vcflib.read(vcf_filename)

    # Figure out the Functional annotations.
    assert vcf.matrix.headerlines, "No header lines"
    x = [x for x in vcf.matrix.headerlines if x.find("<ID=ANN,") >= 0]
    if not x:
        return
    # No duplicates.
    # The ANN line can end with:
    #   ERRORS / WARNINGS / INFO'">
    #   ERRORS / WARNINGS / INFO' ">
    # I encountered a VCF file that contained two ANN lines differing
    # by this spacing.  Normalize these lines and make sure there are
    # no duplicates.
    x = [x.replace("ERRORS / WARNINGS / INFO' \">",
                   "ERRORS / WARNINGS / INFO'\">") for x in x]
    x = {}.fromkeys(x).keys()
    assert len(x) == 1, "Multiple ANN headers: %s" % vcf_filename
    header = x[0]
    x = header.strip()
    TEXT = "Functional annotations:"
    assert TEXT in x
    x = x[x.index(TEXT)+len(TEXT):]   # Get rid of "Functional annotations:"
    assert x.endswith('">')           # No ">
    x = x[:-2].strip()
    assert x.startswith("'") and x.endswith("'")   # No ''
    x = x[1:-1]
    x = x.split("|")
    x = [x.strip() for x in x]
    annotations = x

    handle = open(out_filename, 'w')
    header = ["Chrom", "Pos", "Ref", "Alt"] + annotations
    print >>handle, "\t".join(header)

    for i in range(vcf.num_variants()):
        var = vcflib.get_variant(vcf, i)
        if "ANN" not in var.infodict:
            continue

        # Can have multiple annotations if there are more than one allele.
        # <ALLELE>|...|...|,<ALLELE>|...|...|
        # If this happens, just add them to the file.
        x = var.infodict["ANN"]
        annots = x.split(",")
        for annot in annots:
            x = annot.split("|")
            x = [x.strip() for x in x]
            values = x
            assert len(values) == len(annotations), \
                   "Mismatch annotations %d %d: %s %s %d" % (
                len(annotations), len(values), vcf_filename,
                var.chrom, var.pos)

            alt = ",".join(var.alt)
            x = [var.chrom, var.pos, var.ref, alt] + values
            assert len(x) == len(header)
            print >>handle, "\t".join(map(str, x))


def merge_parsed_files(parsed_files, outfile):
    # First, make sure each of the parsed files has the same header.
    from genomicode import filelib

    assert parsed_files

    header = None
    for f in parsed_files:
        cols = filelib.read_cols(f).next()
        if not header:
            header = cols
        assert header == cols, "Mismatched headers"
    assert header

    handle = open(outfile, 'w')
    seen = {}
    for f in parsed_files:
        for line in filelib.openfh(f):
            if line in seen:
                continue
            seen[line] = 1
            print >>handle, line,


def clean_snpeff_file(infile, outfile):
    import os
    from genomicode import filelib
    from genomicode import sortlib

    if not os.path.exists(infile):
        return

    # Chrom
    # Pos
    # Ref
    # Alt
    # Allele                    Multiple Feature_ID per Allele.
    # Annotation
    # Annotation_Impact
    # Gene_Name                 Unique, given coord and Feature_ID.
    # Gene_ID                   Unique, given coord and Feature_ID.
    # Feature_Type              Unique, given coord and Feature_ID.
    # Feature_ID                Can have multiple Feature_ID per coord.
    # Transcript_BioType        Unique, given coord and Feature_ID.
    # Rank
    # HGVS.c
    # HGVS.p
    # cDNA.pos / cDNA.length
    # CDS.pos / CDS.length
    # AA.pos / AA.length
    # Distance
    # ERRORS / WARNINGS / INFO  Unique, given coord and Feature_ID.

    if False:  # For debugging
        # Check to see which columns are unique given the same coordinate.
        coord2header2values = {}  # coord -> header -> list of values
        for d in filelib.read_row(infile, header=1):
            coord = d.Chrom, d.Pos, d.Ref, d.Alt, d.Feature_ID
            if coord not in coord2header2values:
                coord2header2values[coord] = {}
            header2values = coord2header2values[coord]
            for header in d._nheader:
                if header in ["Chrom", "Pos", "Ref", "Alt", "Allele"]:
                    continue
                value = getattr(d, header)
                if header not in header2values:
                    header2values[header] = []
                if value not in header2values[header]:
                    header2values[header].append(value)
        not_unique = {}  # header -> 1
        for coord, header2values in coord2header2values.iteritems():
            for header, values in header2values.iteritems():
                if len(values) > 1:
                    not_unique[header] = 1
        for x in sorted(not_unique):
            print x

    # Just merge each of these annotations using commas.
    coord2ds = {}  # (chrom, pos, ref, alt) -> list of d
    header = None
    for d in filelib.read_row(infile, header=1):
        if not header:
            header = d._header
        coord = d.Chrom, d.Pos, d.Ref, d.Alt
        coord2ds.setdefault(coord, []).append(d)
    # Sort each record by Allele, Gene_Name, Feature_Type, Feature_ID
    for coord, ds in coord2ds.iteritems():
        schwartz = [
            (d.Allele, d.Gene_Name, d.Feature_Type, d.Feature_ID, d)
            for d in ds]
        schwartz.sort()
        ds = [x[-1] for x in schwartz]
        coord2ds[coord] = ds
    # Convert ds to matrix of values
    coord2matrix = {}  # coord -> list of lists
    for coord, ds in coord2ds.iteritems():
        matrix = [d._cols for d in ds]
        coord2matrix[coord] = matrix
    # Merge the values.
    DELIM = ","
    coord2row = {}    # coord -> list
    for coord, matrix in coord2matrix.iteritems():
        row = []
        # Coordinate should be unique.
        for j in range(4):
            for i in range(len(matrix)):
                assert matrix[i][j] == matrix[0][j]
            row.append(matrix[0][j])
        assert len(matrix)
        for j in range(4, len(matrix[0])):
            x = [matrix[i][j] for i in range(len(matrix))]
            for y in x:
                assert DELIM not in y
            # If every value is the same, then just use the first value.
            z = sorted(x)
            if z[0] == z[-1]:
                x = [z[0]]
            x = ",".join(x)
            row.append(x)
        assert len(row) == len(header)
        coord2row[coord] = row
    # Write out each of the rows.
    h = sortlib.hash_natural
    all_coords = coord2row.keys()
    schwartz = [(h(x[0]), h(x[1]), x[2], x[3], x) for x in all_coords]
    schwartz.sort()
    all_coords = [x[-1] for x in schwartz]
    handle = open(outfile, 'w')
    print >>handle, "\t".join(header)
    for coord in all_coords:
        x = coord2row[coord]
        assert len(x) == len(header)
        print >>handle, "\t".join(map(str, x))


def add_snpeff_to_svm(svm_file, snpeff_file, outfile):
    import shutil
    from genomicode import filelib
    from genomicode import SimpleVariantMatrix
    from genomicode import AnnotationMatrix

    if not filelib.exists_nz(snpeff_file):
        shutil.copy2(svm_file, outfile)
        return

    # Read the annotations.
    header = None  # includes Chrom, Pos, Ref, Alt
    coord2d = {}
    for d in filelib.read_row(snpeff_file, header=1):
        if header is None:
            header = d._header
        coord = d.Chrom, d.Pos, d.Ref, d.Alt
        coord2d[coord] = d
    
    svm = SimpleVariantMatrix.read_as_am(svm_file)
    CHROM = svm.header2annots["______Chrom"]
    POS = svm.header2annots["______Pos"]
    REF = svm.header2annots["______Ref"]
    ALT = svm.header2annots["______Alt"]

    snpeff_header = header[4:]
    snpeff_matrix = []  # Row major.
    for i in range(len(CHROM)):
        coord = CHROM[i], POS[i], REF[i], ALT[i]
        row = [""] * len(snpeff_header)
        d = coord2d.get(coord)
        if d:
            row = d._cols[4:]
        assert len(row) == len(snpeff_header)
        snpeff_matrix.append(row)
    assert len(snpeff_matrix) == len(CHROM)
    # AnnotationMatrix is column major.
    snpeff_annots = []
    for j in range(len(snpeff_header)):
        x = [snpeff_matrix[i][j] for i in range(len(snpeff_matrix))]
        snpeff_annots.append(x)
    # Convert the headers to SVM format.
    snpeff_header = ["SnpEff______%s" % x for x in snpeff_header]
    # Make the new SimpleVariantMatrix.
    headers = svm.headers[:4] + snpeff_header + svm.headers[4:]
    x = [svm.header2annots[x] for x in svm.headers_h]
    all_annots = x[:4] + snpeff_annots + x[4:]
    merged = AnnotationMatrix.create_from_annotations(
        headers, all_annots, headerlines=svm.headerlines)
    SimpleVariantMatrix.write_from_am(outfile, merged)
