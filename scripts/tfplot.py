#!/usr/bin/env python

import sys, os

# Classes:
# BasePair2PixConverter
# SequenceLayout
# TFBSLayout
# LabelLayout
#
# Functions:
# plot
# plot_sequence
# plot_tfbs
#
# resolve_genome_files
# resolve_symbol_or_file
# resolve_sequence
# resolve_matrices
# overlaps
# sign
# safe_split_int

class BasePair2PixConverter:
    # Calculates coordinates from base pairs to pixels.
    # - Specify origin of base pairs (x_bp, y_bp).
    # - Specify corresponding origin of pixels (x_pix, y_pix).
    # - Y axis of pixels starts at top and increases downwards.
    # - Y axis of base pairs starts at bottom and increases upwards.
    # - X axis of pixels starts at left and increases rightwards.
    # - X axis of base pairs increases or decreases depending on the
    #   strand of the transcript.
    def __init__(self, pixels_per_bp, x_bp, y_bp, x_pix, y_pix, strand):
        self.pixels_per_bp = pixels_per_bp
        self.xoff_bp, self.yoff_bp = x_bp, y_bp
        self.xoff_pix, self.yoff_pix = x_pix, y_pix
        self.strand = strand
    def pixels(self, x_bp, y_bp):
        # Return coordinate in pixels.
        import math
        if self.strand == "+":
            x_bp = x_bp - self.xoff_bp
        else:
            x_bp = self.xoff_bp - x_bp
        y_bp = y_bp - self.yoff_bp
        x = int(math.ceil(x_bp * self.pixels_per_bp))
        y = int(math.ceil(y_bp * self.pixels_per_bp))
        x_pix = self.xoff_pix + x
        y_pix = self.yoff_pix - y
        return x_pix, y_pix

class SequenceLayout:
    def __init__(self, bp2pix, line_width, tss_height, tss_width,
                 hash_interval, hash_height, color):
        self.bp2pix = bp2pix
        self.line_width = line_width
        self.tss_height = tss_height
        self.tss_width = tss_width
        self.hash_interval = hash_interval
        self.hash_height = hash_height
        self._color = color
    def vthicken(self, x, y, width, height, num_pixels):
        # Thicken a vertical line up to a specific number of pixels.
        import math
        np = num_pixels - width
        if np <= 0:
            return x, y, width, height
        hnp = int(math.floor(np/2.0))
        return x-hnp, y, width+np, height
    def hthicken(self, x, y, width, height, num_pixels):
        # Thicken a horizontal line up to a specific number of pixels.
        import math
        np = num_pixels - height
        if np <= 0:
            return x, y, width, height
        hnp = int(math.floor(np/2.0))
        return x, y-hnp, width, height+np
    def sequence_line(self, start, length):
        # Get the x, y, width, height coordinates of a sequence.
        x1, y1 = self.bp2pix.pixels(start, 0)
        x2, y2 = self.bp2pix.pixels(start+length, 0)
        width, height = x2-x1, y2-y1
        x = self.hthicken(x1, y1, width, height, self.line_width)
        return x
    def sequence_hash_lines(self, start, length, tss):
        # Return a list of (x, y, width, height) to plot hash marks
        # along the sequence.
        # Adjust the start and length to make sure each hash mark is
        # covered.
        if self.bp2pix.strand == "+":
            start = start - self.hash_interval
        else:
            start = start + self.hash_interval
        length = length + self.hash_interval*2

        # Figure out where to draw the hash marks.
        hash_start = start/self.hash_interval*self.hash_interval
        if tss is not None:
            # If the TSS is specified, draw relative to the TSS.
            hash_start += tss % self.hash_interval
        hash_end = start + length
            
        lines = []
        # Round to the nearest hash interval
        for pos in range(hash_start, hash_end, self.hash_interval):
            x1, y1 = self.bp2pix.pixels(pos, -self.hash_height/2)
            x2, y2 = self.bp2pix.pixels(pos, self.hash_height/2)
            width, height = x2-x1, y2-y1
            x = x1, y1, width, height
            lines.append(x)
        return lines
    
    def tss_coords(self, tss):
        # Return an array of 5 coordinates (x, y):
        # 1.  Bottom of ascending line.
        # 2.  Top of ascending line.
        # 3.  Point of arrow.
        # 4.  Top tail of arrow.
        # 5.  Bottom tail of arrow.
        if not tss:
            return []
        
        ARROW_BP = self.tss_height*0.375
        x1, y1 = self.bp2pix.pixels(tss, 0)
        x2, y2 = self.bp2pix.pixels(tss, self.tss_height)
        if self.bp2pix.strand == "+":
            x3, y3 = self.bp2pix.pixels(tss+self.tss_width, self.tss_height)
        else:
            x3, y3 = self.bp2pix.pixels(tss-self.tss_width, self.tss_height)
        if self.bp2pix.strand == "+":
            x4, y4 = self.bp2pix.pixels(
                tss+self.tss_width-ARROW_BP, self.tss_height+ARROW_BP)
            x5, y5 = self.bp2pix.pixels(
                tss+self.tss_width-ARROW_BP, self.tss_height-ARROW_BP)
        else:
            x4, y4 = self.bp2pix.pixels(
                tss-self.tss_width+ARROW_BP, self.tss_height+ARROW_BP)
            x5, y5 = self.bp2pix.pixels(
                tss-self.tss_width+ARROW_BP, self.tss_height-ARROW_BP)
        coords = [(x1, y1), (x2, y2), (x3, y3), (x4, y4), (x5, y5)]
        return coords
    def color(self):
        return self._color

class TFBSLayout:
    def __init__(self, bp2pix, tfbs_height, matrix2color, matrix2style,
                 color_gain):
        self.bp2pix = bp2pix
        self.tfbs_height = tfbs_height
        self.matrix2color = matrix2color
        self.matrix2style = matrix2style
        self.color_gain = color_gain
    def coord(self, matrix, pos, length, strand, NLP, height_offset):
        x1_bp, x2_bp = pos, pos+length
        y1_bp, y2_bp = 0, self.tfbs_height
        if strand != self.bp2pix.strand:
            y2_bp = -self.tfbs_height
        x1, y1 = self.bp2pix.pixels(x1_bp, y1_bp)
        x2, y2 = self.bp2pix.pixels(x2_bp, y2_bp)
        if strand != self.bp2pix.strand:
            y1, y2 = y1+1, y2+1
        width, height = x2-x1, y2-y1
        height += height_offset*sign(height)
        return x1, y1, width, height
    def _scale_color(self, color, percent):
        assert percent >= 0 and percent <= 1
        # Apply a default scaling factor to the color.
        DEF_SCALE_PERCENT = 2
        percent = min(percent*DEF_SCALE_PERCENT*self.color_gain, 1.0)
        
        r, g, b = color
        r = int(round(255-(255-r)*percent))
        g = int(round(255-(255-g)*percent))
        b = int(round(255-(255-b)*percent))
        #r = int(round(r*percent))
        #g = int(round(g*percent))
        #b = int(round(b*percent))
        return r, g, b
    def color(self, matrix, pos, length, strand, NLP):
        import math
        DEFAULT = 0, 0, 0

        # NLP is log_e(p-value).
        min_NLP = -math.log(0.05)
        max_NLP = -math.log(0.0000001)
        percent = (NLP-min_NLP) / (max_NLP-min_NLP)
        percent = max(min(percent, 1.0), 0.0)

        color = self.matrix2color.get(matrix, DEFAULT)
        color = self._scale_color(color, percent)
        return color
    def style(self, matrix, pos, length, strand, NLP):
        DEFAULT = "f"
        if matrix not in self.matrix2style:
            return DEFAULT
        style = self.matrix2style[matrix]
        assert style in "fo"
        return style
    
class LabelLayout:
    def __init__(self, bp2pix, label_height):
        self.bp2pix = bp2pix
        self.height = label_height
        self._color = 0, 0, 0
    def coord(self):
        OFFSET = 3
        if self.bp2pix.strand == "-":
            OFFSET = -OFFSET
        x, y = self.bp2pix.pixels(self.bp2pix.xoff_bp+OFFSET, self.height)
        return x, y
    def color(self):
        return self._color

def plot(plotlib, image, seq_layout, tfbs_layout, label_layout, 
         gene_symbol, gen_start, gen_length, tss, binding_sites):
    # If TFBS overlap, then increase the height of one of them so that
    # you can still distinguish independent binding sites.
    # Sort them by strand, position.
    # If the sequence is on the "-" strand, need to reverse the
    # sorting of the position.
    s = 1
    if seq_layout.bp2pix.strand == "-":
        s = -1
    schwartz = sorted([(x[3], x[1]*s, x) for x in binding_sites])
    binding_sites = [x[-1] for x in schwartz]
    offsets = [0] * len(binding_sites)
    for i in range(1, len(binding_sites)):
        bs1 = binding_sites[i-1]
        bs2 = binding_sites[i]
        pos1, len1 = bs1[1], bs1[2]
        pos2, len2 = bs2[1], bs2[2]
        s1, e1 = pos1, pos1+len1
        s2, e2 = pos2, pos2+len2

        if s2 >= s1 and e2 <= e1:
            # If this TFBS is contained within the previous one, then
            # make them the same height.
            offsets[i] = offsets[i-1]
        elif overlaps((s1, e1), (s2, e2)):
            offsets[i] = offsets[i-1]+1
    
    # Draw each of the transcription factor binding sites.  Sort by
    # increasing NLP so that more significant ones show up front.
    schwartz = sorted([(x[4], x, y) for (x, y) in zip(binding_sites, offsets)])
    binding_sites = [x[-2] for x in schwartz]
    offsets = [x[-1] for x in schwartz]

    # Draw each of the binding sites.
    for i in range(len(binding_sites)):
        tfbs = binding_sites[i]
        offset = offsets[i]
        matrix, pos, length, strand, nlp = tfbs
        plot_tfbs(
            plotlib, image, tfbs_layout, matrix, pos, length, strand, nlp,
            offset)
        
    # Draw the underlying sequence.
    plot_sequence(plotlib, image, seq_layout, gen_start, gen_length, tss)

    if label_layout is not None:
        plot_label(plotlib, image, label_layout, gene_symbol)

def plot_sequence(plotlib, image, seq_layout, start, length, tss):
    COLOR = seq_layout.color()

    # Draw the sequence.
    x = seq_layout.sequence_line(start, length)
    x, y, width, height = x
    plotlib.line(image, x, y, width, height, COLOR)

    # Draw the TSS.
    #for x in seq_layout.tss_lines(tss):
    #    x, y, width, height = x
    #    plotlib.line(image, x, y, width, height, COLOR)
    coords = seq_layout.tss_coords(tss)
    
    # Draw the ascending line.
    (x1, y1), (x2, y2), (x3, y3) = coords[0], coords[1], coords[2]
    width, height = x2-x1, y2-y1
    # The coordinates should have exclusive ends, so add 1 to the
    # width and heights.
    plotlib.line(image, x1, y1, width+sign(width), height-sign(height), COLOR)

    # Draw the horizontal line.
    width, height = x3-x2, y3-y2
    plotlib.line(image, x2, y2, width+sign(width), height+sign(height), COLOR)

    # Draw the arrow head.
    plotlib.polygon(image, coords[2:], COLOR)

    # Draw the hash lines on the sequence.
    for x in seq_layout.sequence_hash_lines(start, length, tss):
        x, y, width, height = x
        plotlib.line(
            image, x, y, width+sign(width), height+sign(height), COLOR)
        
def plot_tfbs(plotlib, image, tfbs_layout, matrix, pos, length, strand, nlp,
              height_offset):
    COLOR = tfbs_layout.color(matrix, pos, length, strand, nlp)
    STYLE = tfbs_layout.style(matrix, pos, length, strand, nlp)
    x = tfbs_layout.coord(matrix, pos, length, strand, nlp, height_offset)
    x, y, width, height = x

    color, outline = COLOR, None
    if STYLE == "o":
        color, outline = None, COLOR
    plotlib.rectangle(image, x, y, width, height, color=color, outline=outline)

def plot_label(plotlib, image, label_layout, gene_symbol):
    COLOR = label_layout.color()
    x, y = label_layout.coord()
    height = int(round(label_layout.height*label_layout.bp2pix.pixels_per_bp))
    fontsize = plotlib.fit_fontsize_to_height(height)
    plotlib.text(image, x, y, gene_symbol, fontsize, COLOR)


def resolve_genome_files(genome):
    from genomicode import config

    ra_chrom_var = "RA_CHROM_%s" % genome.upper()
    knowngene_var = "knowngene_%s" % genome.lower()

    errors = []
    if not hasattr(config, ra_chrom_var):
        x = "The ra_chrom file for %s is not configured." % genome
        errors.append(x)
    if not hasattr(config, knowngene_var):
        x = "The knowngene file for %s is not configured." % genome
        errors.append(x)
    assert not errors, "\n".join(errors)

    ra_chrom_path = getattr(config, ra_chrom_var)
    knowngene_file = getattr(config, knowngene_var)
    if not os.path.exists(ra_chrom_path):
        x = "I could not find the ra_chrom path: %s." % ra_chrom_path
        errors.append(x)
    if not os.path.exists(knowngene_file):
        x = "I could not find the knowngene file: %s." % knowngene_file
        errors.append(x)
    assert not errors, "\n".join(errors)

    x = ra_chrom_path, knowngene_file
    return x


def resolve_symbol_or_file(name):
    from genomicode import filelib
    
    if not os.path.exists(name):
        return [name]
    symbols = [x.strip() for x in filelib.openfh(name)]
    return symbols

def resolve_sequence(
    name, bp_upstream, length, ra_chrom_path, knowngene_file, 
    default_transcript=None, skip_unknown_genes=False):
    from genomicode import parselib
    from genomicode import genomelib

    transcript_num = default_transcript
    gene_symbol = name
    if "," in name:
        gene_symbol, transcript_num = name.split(",", 1)
        assert transcript_num, "Invalid gene symbol: %r" % name
        transcript_num = int(transcript_num)

    proms = genomelib.get_promoters(
        gene_symbol, -bp_upstream, length, gene_file=knowngene_file,
        ra_path=ra_chrom_path)
    if not proms and skip_unknown_genes:
        return None
    assert proms, "I could not find gene: %s" % gene_symbol
    
    # If there is only 1 promoter, then use that one.
    if len(proms) == 1 and transcript_num is None:
        transcript_num = 0

    if len(proms) > 1 and transcript_num is None:
        print "Multiple transcripts for %s." % gene_symbol
        print "Please specify a specific transcript."
        x = "Index", "Chrom", "Strand", "TSS"
        print "\t".join(map(str, x))
        for i in range(len(proms)):
            chrom, tss, strand, start, length, seq = proms[i]
            tss = parselib.pretty_int(tss)
            x = i, chrom, strand, tss
            print "\t".join(map(str, x))
        sys.exit(0)

    assert transcript_num is not None
    assert transcript_num >= 0 and transcript_num < len(proms), \
           "Invalid transcript: %s" % transcript_num
    
    x = proms[transcript_num]
    chrom, tss, txn_strand, prom_base, prom_length, prom_seq = x
    return gene_symbol, chrom, prom_base, prom_length, txn_strand, tss

def resolve_matrices(names, all_matrices):
    # Clean up user input.
    from genomicode import motiflib

    if all_matrices:
        names = [x[0] for x in motiflib.list_matrices()]

    # Parse out the matrices from the user input.
    matrices = []
    matrix2color = {}
    matrix2style = {}
    for name in names:
        # Split the name from the color.
        color = style = None
        x = name.split(",")
        if len(x) == 2:
            name, color = x
        elif len(x) == 3:
            name, color, style = x
        assert len(x) < 4

        # Map name to matrix IDs.
        if motiflib.is_matrix_id(name):
            matrix_ids = [name]
        else:
            x = motiflib.gene2matrices(name)
            matrix_ids = [x.matid for x in x]
        assert matrix_ids, "Unrecognized gene or matrix name: %s" % name
        matrices.extend(matrix_ids)

        # Map the matrix IDs to the color.
        if color is not None:
            assert color.startswith("0x"), "invalid color format: %s" % color
            assert len(color) == 8, "invalid color format: %s" % color
            r, g, b = color[2:4], color[4:6], color[6:8]
            r, g, b = int(r, 16), int(g, 16), int(b, 16)
            for matid in matrix_ids:
                matrix2color[matid] = r, g, b

        if style is not None:
            assert style in "fo"
            for matid in matrix_ids:
                matrix2style[matid] = style
        
    matrices = sorted({}.fromkeys(matrices))
    return matrices, matrix2color, matrix2style

def overlaps(r1, r2):
    if r1 > r2:
        r1, r2 = r2, r1
    s1, e1 = r1
    s2, e2 = r2
    assert s1 <= e1
    assert s2 <= e2
    return e1 > s2

def sign(x):
    if x < 0:
        return -1
    return 1

def safe_split_int(s):
    return [int(x) for x in s.split(",") if x]

def main():
    from optparse import OptionParser, OptionGroup

    import math
    from genomicode import genomelib
    from genomicode import motiflib
    from genomicode import parselib
    
    usage = "usage: %prog [options] <GENE SYMBOL>[,#] [...]"
    parser = OptionParser(usage=usage, version="%prog 01")

    # Matrix options.
    parser.add_option(
        "-m", "--matrix", dest="matrices", default=[], action="append",
        help="Plot binding sites for this matrix or gene.  "
        "Format: <matrix>[,<color>[,<f|o>].  Color is 0xRRGGBB format.  "
        "f(illed)|o(outline); def=f.")
    parser.add_option(
        "--all_matrices", default=False, action="store_true",
        help="Search all matrices.")

    # Gene options.
    parser.add_option(
        "--genome", default="hg19",
        help="Genome to search, e.g. hg19 (default), mm11.")
    parser.add_option(
        "--upstream", dest="upstream", default=250, type="int",
        help="Number of base pairs upstream of the TSS (def 250).")
    parser.add_option(
        "--downstream", dest="downstream", default=50, type="int",
        help="Number of base pairs downstream of the TSS (def 50).")
    parser.add_option(
        "-s", dest="strict", default=False, action="store_true",
        help="Use strict checking of gene names.")

    # Plotting options.
    parser.add_option(
        "--output_as_table", default=False, action="store_true")
    parser.add_option(
        "--pvalue", dest="pvalue_cutoff", default=None, type="float",
        help="p-value cutoff.")
    parser.add_option(
        "--max_genes", dest="max_genes", default=None, type="int",
        help="Maximum number of genes to plot.")
    parser.add_option(
        "-l", dest="label", default=False, action="store_true",
        help="Label the gene names.")
    parser.add_option(
        "-g", dest="gain", default=1, type="float",
        help="Increase the gain of the colors of the TFBS (def 1).")
    parser.add_option(
        "--format", dest="image_format", type="choice",
        choices=["png", "svg"], default="png",
        help="Image format: png (default) or svg.")

    # Running options.
    parser.add_option(
        "-j", "--jobname", dest="jobname", type="string", default="out",
        help="Name of output file.")
    parser.add_option(
        "--num_procs", dest="num_procs", type="int", default=1,
        help="Number of jobs to run in parallel.")

    options, args = parser.parse_args()
    if not args:
        parser.error("Please specify a gene to analyze.")
    if options.num_procs < 1 or options.num_procs > 100:
        parser.error("Please specify between 1 and 100 processes.")

    ra_chrom_path, knowngene_file = resolve_genome_files(options.genome)
         
    gene_symbols = []
    for symbol_or_file in args:
        x = resolve_symbol_or_file(symbol_or_file)
        gene_symbols.extend(x)

    default_transcript = None
    skip_unknown_genes = False
    if not options.strict:
        default_transcript = 0
        skip_unknown_genes = True

    if not options.matrices and not options.all_matrices:
        parser.error("Please specify a matrix to plot.")
    if options.upstream < 0:
        parser.error("Upstream should be 0 or positive.")
    if options.downstream < 0:
        parser.error("Downstream should be 0 or positive.")
    seq_length = options.upstream + options.downstream
    if seq_length <= 0:
        parser.error("No sequence.")

    # Choose a plotting library.
    if options.jobname.lower().endswith(".png"):
        options.jobname = options.jobname[:-4]
    if options.image_format == "svg":
        plotlib = __import__(
            "genomicode.svgplot", globals(), locals(), ["svgplot"])
        outfile = "%s.svg" % options.jobname
    else:
        plotlib = __import__(
            "genomicode.pilplot", globals(), locals(), ["pilplot"])
        outfile = "%s.png" % options.jobname

    # Figure out the genes to plot.
    GENES = []   # list of gene_symbol, chrom, start, length, strand, tss
    for name in gene_symbols:
        # Figure out where the gene lies on the chromosome.
        x = resolve_sequence(
            name, options.upstream, seq_length, ra_chrom_path, knowngene_file,
            default_transcript=default_transcript,
            skip_unknown_genes=skip_unknown_genes)
        if x is None:
            continue
        GENES.append(x)

    if options.max_genes is not None:
        if options.max_genes <= 0:
            parser.error("Invalid max_genes.")
        GENES = GENES[:options.max_genes]

    # Figure out the matrices to search for.
    x = resolve_matrices(options.matrices, options.all_matrices)
    matrices, matrix2color, matrix2style = x

    # Constants governing how the figures are drawn.
    PIXELS_PER_BP = 4         # How many pixels for each base pair.
    X_BORDER = 5              # Number of pixels of border around the figure.
    Y_BORDER = 5
    GENE_BUFFER = 4           # Buffer between genes, in pixels.
    
    SEQ_COLOR = (0, 0, 0)     # Color for the sequence.
    SEQ_WIDTH = 1             # Width of line for sequence, in pixels.
    TSS_HEIGHT = 3            # Height of the TSS, in base pairs.
    TSS_WIDTH = TSS_HEIGHT*2  # Width of the TSS
    TFBS_HEIGHT = 2           # Height of TFBS, in base pairs.
    HASH_INTERVAL = 50        # Number of base pairs per hash mark.
    HASH_HEIGHT = 0.5         # Height of hash marks, in base pairs
    LABEL_HEIGHT = 5          # Height of label, in base pairs.

    NUM_GENES = len(GENES)

    # Figure out the geometry for a single gene.
    plot_width = int(PIXELS_PER_BP * seq_length)
    # Hack: Add some to TFBS_HEIGHT to handle overlapping TFBS.
    plot_height = int(PIXELS_PER_BP * max(TSS_HEIGHT, TFBS_HEIGHT+6)*2)
    # Figure out the geometry for the entire figure.
    total_width = X_BORDER*2 + plot_width
    total_height = (
        Y_BORDER*2 + plot_height*NUM_GENES + GENE_BUFFER*(NUM_GENES-1))

    image = plotlib.image(total_width, total_height)

    if options.output_as_table:
        header = [
            "Gene Symbol", "Chromosome", "Strand", "Transcription Start",
            "Matrix ID", "Matrix Gene Symbol", "TFBS Pos", "TFBS Strand",
            "TFBS Offset", "P-value", "Sequence"]
        print "\t".join(header)

    # Plot each gene.
    for gene_num, x in enumerate(GENES):
        gene_symbol, chrom, gen_start, gen_length, txn_strand, tss = x
        if not options.output_as_table:
            print "%s: %+d to %+d relative to TSS at chr%s:%s:%s." % (
                gene_symbol, -options.upstream, -options.upstream+seq_length,
                chrom, parselib.pretty_int(tss), txn_strand)

        # Get the TFBS from that site.
        nlp_cutoff = 0
        if options.pvalue_cutoff:
            nlp_cutoff = -math.log(options.pvalue_cutoff)
        data = motiflib.score_tfbs_genome(
            chrom, gen_start, gen_length, matrices=matrices, nlp=nlp_cutoff,
            num_procs=options.num_procs, ra_path=ra_chrom_path)
        
        ## # If multiple matrices for the same gene symbol hit the same
        ## # place, then keep the one with the highest NLP.
        ## # Sort by gene_symbol, position, decreasing NLP.
        ## x = [(matid2info[x[0]].Gene_Symbol, x[3], -x[4], x) for x in data]
        ## x.sort()
        ## data = [x[-1] for x in x]
        ## i = 0
        ## while i < len(data)-1:
        ##     di, dj = data[i], data[i+1]
        ##     gs0 = matid2info[di[0]].Gene_Symbol
        ##     gs1 = matid2info[dj[0]].Gene_Symbol
        ##     # Same matrix and position.
        ##     if gs0 == gs1 and di[3] == dj[3]:
        ##         del data[i+1]
        ##     else:
        ##         i += 1

        # Sort by position relative to TSS, decreasing NLP
        x = [(genomelib.genbase2tssbase(x[3], tss, txn_strand), -x[4], x)
             for x in data]
        x.sort()
        data = [x[-1] for x in x]

        # Inspect the data to improve the formatting of the results.
        name_len = 0
        gs_len = 0
        tss_dist_len = 0
        #NLP_len = 0
        for x in data:
            matrix, chrom, strand, pos, NLP = x

            m = motiflib.matid2matrix(matrix)
            
            name_len = max(name_len, len(matrix))
            gs_len = max(gs_len, len(m.gene_symbol))
            tss_dist = genomelib.calc_tss_seq_dist(
                tss, txn_strand, pos, m.length)
            tss_dist_len = max(
                tss_dist_len, len(parselib.pretty_int(tss_dist)))
            #NLP_len = max(NLP_len, len("%.2f" % NLP))

        # Print the data.
        binding_sites = []
        for x in data:
            matrix, chrom, strand, pos, NLP = x
            
            m = motiflib.matid2matrix(matrix)
            
            x = matrix, pos, m.length, strand, NLP
            binding_sites.append(x)

            # Get the matrix sequence, with some flanking bases.
            FLANK = int(math.ceil((20-m.length)/2.0))
            FLANK = max(FLANK, 2)
            left_flank = min(FLANK, pos)
            right_flank = FLANK
            seq_pos = pos - left_flank
            seq_len = m.length + left_flank + right_flank
            seq = genomelib.get_sequence(
                chrom, seq_pos, seq_len, ra_path=ra_chrom_path)
            s1 = seq[:left_flank]
            s2 = seq[left_flank:left_flank+m.length]
            s3 = seq[left_flank+m.length:]
            seq = s1.lower() + s2.upper() + s3.lower()
            
            # Print to stdout.
            mat_strand = "+"
            if strand != txn_strand:
                mat_strand = "-"
            tss_dist = genomelib.calc_tss_seq_dist(
                tss, txn_strand, pos, m.length)
            pvalue = parselib.pretty_pvalue(math.exp(-NLP), nsig=2)
            if options.output_as_table:
                x = (gene_symbol, chrom, txn_strand, tss,
                     matrix, m.gene_symbol,
                     pos, strand, tss_dist, pvalue, seq)
                assert len(x) == len(header)
                print "\t".join(map(str, x))
            else:
                #print "  %-*s [%*s] %s:%s:%s (%*s:%s) %*.2f %s" % (
                print "  %-*s [%*s] %s:%s (%*s) %-8s %s" % (
                    gs_len, m.gene_symbol, name_len, matrix, 
                    parselib.pretty_int(pos), strand, tss_dist_len,
                    parselib.pretty_int(tss_dist), pvalue, seq)

        # Initialize some objects to plot the figures.
        bp_start = gen_start
        if txn_strand == "-":
            bp_start = gen_start+gen_length
        x_bp = bp_start
        y_bp = 0
        x_pix = X_BORDER
        y_pix = Y_BORDER + gene_num*(plot_height+GENE_BUFFER) + plot_height/2
    
        bp2pix = BasePair2PixConverter(
            PIXELS_PER_BP, x_bp, y_bp, x_pix, y_pix, txn_strand)
        seq_layout = SequenceLayout(
            bp2pix, SEQ_WIDTH, TSS_HEIGHT, TSS_WIDTH,
            HASH_INTERVAL, HASH_HEIGHT, SEQ_COLOR)
        tfbs_layout = TFBSLayout(
            bp2pix, TFBS_HEIGHT, matrix2color, matrix2style, options.gain)
        label_layout = None
        if options.label:
            label_layout = LabelLayout(bp2pix, LABEL_HEIGHT)

        # Actually plot the figures.
        plot(
            plotlib, image, seq_layout, tfbs_layout, label_layout,
            gene_symbol, gen_start, gen_length, tss, binding_sites)
        sys.stdout.flush()

    plotlib.write(image, open(outfile, 'w'))
    
if __name__ == '__main__':
    main()
