from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_filename):
        #import shutil
        from genomicode import filelib
        from genomicode import parallel
        from genomicode import alignlib
        from genomicode import SimpleVariantMatrix
        from genomicode import AnnotationMatrix
        from Betsy import module_utils as mlib

        summary_node = in_data
        summary_filename = summary_node.identifier
        metadata = {}

        buildver = mlib.get_user_option(
            user_options, "annovar_buildver", allowed_values=["hg19"],
            not_empty=True)

        # Name files.
        p, root, ext = mlib.splitpath(summary_filename)
        annovar_infile = "pos.txt"
        log_filename = "annovar.log"
        # Annovar takes a filestem, without the ".vcf".
        annovar_outstem = "annotations"
        # Produces file:
        # <annovar_outstem>.hg19_multianno.txt
        multianno_file = "%s.hg19_multianno.txt" % annovar_outstem
        #temp_file = "temp.txt"

        # Make the infile for Annovar.
        # <chrom> <start> <end> <ref> <alt>
        handle = open(annovar_infile, 'w')
        for d in filelib.read_row(summary_filename, skip=2, header=1):
            x = d.Chrom, d.Pos, d.Pos, d.Ref, d.Alt
            print >>handle, "\t".join(x)
        handle.close()

        cmd = alignlib.make_annovar_command(
            annovar_infile, log_filename, annovar_outstem, buildver,
            vcf_input=False)
        parallel.sshell(cmd)
        metadata["commands"] = [cmd]

        filelib.assert_exists_nz(log_filename)
        filelib.assert_exists_nz(multianno_file)

        matrix = SimpleVariantMatrix.read(summary_filename)
        annot_matrix = matrix.annot_matrix
        #headers = annot_matrix.headers + anno_header[5:]
        chrom, pos = annot_matrix["Chrom"], annot_matrix["Pos"]
        ref, alt = annot_matrix["Ref"], annot_matrix["Alt"]
        pos = [int(x) for x in pos]

        # Read in the multianno output file.
        pos2d = {}  # (chrom, start, ref, alt) -> d
        anno_header = None
        for d in filelib.read_row(multianno_file, header=1):
            key = d.Chr, int(d.Start), d.Ref, d.Alt
            assert key not in pos2d, "Duplicate pos: %s" % str(key)
            pos2d[key] = d
            if not anno_header:
                anno_header = d._header
        assert anno_header

        # Multianno starts with:
        # Chr Start End Ref Alt
        # Ignore these.
        assert anno_header[:5] == ["Chr", "Start", "End", "Ref", "Alt"]
        headers = anno_header[5:]
        
        all_annots = []
        #for h in annot_matrix.headers_h:
        #    x = annot_matrix.header2annots[h]
        #    all_annots.append(x)
        for i in range(5, len(anno_header)):
            annots = []
            for coord in zip(chrom, pos, ref, alt):
                d = pos2d.get(coord)
                x = ""
                if d:
                    x = d._cols[i]
                annots.append(x)
            all_annots.append(annots)
        x = AnnotationMatrix.create_from_annotations(headers, all_annots)
        matrix.named_matrices.append(("Annovar", x))
        
        SimpleVariantMatrix.write(out_filename, matrix)
        
        ## cols_to_add = len(anno_header) - 5
        ## assert cols_to_add > 0

        ## # Merge the multianno file with the simple call summary.  Add
        ## # these columns before the <Sample>.
        ## # Sample                <Sample>
        ## # Caller                <Caller>
        ## # Chrom  Pos  Ref  Alt  Ref/Alt/VAF
        ## handle = open(temp_file, 'w')
        ## it = filelib.read_cols(summary_filename)
        ## header1 = it.next()
        ## header2 = it.next()
        ## header3 = it.next()
        ## assert len(header1) == len(header2), "%d %d %d %s" % (
        ##     len(header1), len(header2), len(header3), summary_filename)
        ## assert len(header1) == len(header3), "%d %d %d %s" % (
        ##     len(header1), len(header2), len(header3), summary_filename)
        ## assert header1[0] == "Sample"
        ## assert header2[0] == "Caller"
        ## assert header3[:4] == ["Chrom", "Pos", "Ref", "Alt"]
        ## header1 = header1[:4] + [""]*cols_to_add + header1[4:]
        ## header2 = header2[:4] + [""]*cols_to_add + header2[4:]
        ## header3 = header3[:4] + anno_header[5:] + header3[4:]
        ## print >>handle, "\t".join(header1)
        ## print >>handle, "\t".join(header2)
        ## print >>handle, "\t".join(header3)
        ## for cols in it:
        ##     chrom, pos, ref, alt = cols[:4]
        ##     pos = int(pos)
        ##     d = pos2d.get((chrom, pos))
        ##     if not d:
        ##         cols = cols[:4] + [""]*cols_to_add + cols[4:]
        ##         continue
        ##     assert ref == d.Ref, "%s %s %s %s %s %s" % (
        ##         chrom, pos, ref, alt, d.Ref, d.Alt)
        ##     assert alt == d.Alt, "%s %s %s %s %s %s" % (
        ##         chrom, pos, ref, alt, d.Ref, d.Alt)
        ##     x = d._cols[5:]
        ##     assert len(x) == cols_to_add
        ##     cols = cols[:4] + x + cols[4:]
        ##     print >>handle, "\t".join(cols)
        ## handle.close()

        ## shutil.move(temp_file, out_filename)

        return metadata

    
    def name_outfile(self, antecedents, user_options):
        return "calls.annotated.txt"
