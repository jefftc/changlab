"""

Functions:
blat
in_silico_pcr

parse_blat_psl
calc_score_psl

pslCalcMilliBad     Percent identity is: 100.0 - pslCalcMilliBad(psl, True)*0.1
pslScore            Score for PSL output.
pslIsProtein

"""
# External dependencies:
# blat  (binary)


def blat(sequence, delay=60):
    # sequence can either be a single sequence, or a FASTA-formatted
    # list of multiple sequences.  Up to 25 sequences can be submitted
    # at the same time.
    # 
    # Return list of hits as objects with members:
    #  1.  match          Number of base pairs that match that aren't repeats.
    #  2.  mis__match     Base pairs that don't match.
    #  3.  rep__match     Matches that are repeats.
    #  4.  N_s            "N" bases.
    #  5.  Q_gap_count    Number of inserts in query.
    #  6.  Q_gap_bases    Number of bases inserted in query.
    #  7.  T_gap_count
    #  8.  T_gap_bases
    #  9.  strand         Strand of genome.
    # 10.  Q_name         Your name for the sequence.
    # 11.  Q_size
    # 12.  Q_start
    # 13.  Q_end
    # 14.  T_name         Name of target.
    # 15.  T_size
    # 16.  T_start        Start of hit (0-based).
    # 17.  T_end          End of hit (exclusive).
    # 18.  block_count    Blocks (ungapped sequence) in alignment.
    # 19.  blockSizes     comma-separated list
    # 20.  qStarts        comma-separated list
    # 21.  tStarts        comma-separated list
    # o Each variable is a string.
    # o Q_start and Q_end are always relative to the + strand.
    # o qStarts coordinates are reversed on the - strand.  To convert
    #   to + strand coordinates:
    #   qStart = qSize - revQEnd
    #   qEnd = qSize - revQStart
    # o tStarts seems always to be relative to the + strand.
    
    from StringIO import StringIO
    import urllib
    import parselib
    import timer

    timer.wait(delay, name="UCSC")

    URL = "http://genome.ucsc.edu/cgi-bin/hgBlat"

    # For psl format, returns 0-based, exclusive coordinates.  For
    # hyperlink format, returns 1-based inclusive coordinates.
    params = {
        "org" : "Human",
        "db" : "hg18",
        "type" : "BLAT's guess",
        "sort" : "query,score",
        "output" : "psl",
        "userSeq" : sequence,
        # Parameters that appear to be optional.
        #"hgsid" : "117849797",
        #"changInfo" : "",
        #"Submit" : "submit",
        }
    params_str = urllib.urlencode(params)
    url = "%s?%s" % (URL, params_str)
    handle = urllib.urlopen(url)
    data = handle.read()

    x = parselib.get_first_contents(data, "pre")
    if x is None:
        assert data.find("No matches") >= 0, data
        return
    x = parselib.remove_all_tags(x).strip()
    handle = StringIO(x)
    ds = [d for d in parse_blat_psl(handle)]
    return ds


def in_silico_pcr(forward_primer, reverse_primer, delay=5):
    # yield chrom, pos1, pos2
    # pos1 and pos2 define the span of the chromosome that would be
    # amplified in a PCR.
    from StringIO import StringIO
    import urllib
    import parselib
    import timer

    timer.wait(delay, name="UCSC")
    
    URL = "http://genome.ucsc.edu/cgi-bin/hgPcr"

    params = {
        "org" : "Human",
        "db" : "hg18",
        "wp_target" : "genome",
        "wp_f" : forward_primer,
        "wp_r" : reverse_primer,
        "wp_size" : 4000,   # product size
        "wp_perfect" : 15,  # Min Perfect Match
        "wp_good" : 15,     # Min Good Match
        }
    params_str = urllib.urlencode(params)
    url = "%s?%s" % (URL, params_str)
    handle = urllib.urlopen(url)
    data = handle.read()

    x = parselib.get_first_contents(data, "pre")
    if x is None:
        assert data.find("No matches") >= 0, data
        return
    x = parselib.remove_all_tags(x)

    for line in StringIO(x):
        if not line.startswith(">"):
            continue
        x = line.split()
        chr_data = x[1]

        # chr12:61435794+6143584
        # chr11:87781454-87781503
        # "-" means that the PCR product is on the - strand.
        i = chr_data.index(":")
        chrom = chr_data[:i]
        location = chr_data[i+1:]
        i = location.find("+")
        if i < 0:
            i = location.index("-")
        pos1, pos2 = int(location[:i]), int(location[i+1:])
        assert pos1 < pos2
        yield chrom, pos1, pos2

def _clean_blat_psl(handle):
    # Yields cleaned up lines.
    import filelib
    
    handle = filelib.openfh(handle)

    # Format:
    # psLayout version 4 DNA DNA   (OR psLayout version 3)
    #
    # <header 1>
    # <header 2>
    # ----------
    # hits
    assert handle.readline().startswith("psLayout version")
    assert handle.readline().strip() == ""
    # Read the 2 header lines and join them together.
    # header2 contains fewer columns than header1.
    header1 = handle.readline().rstrip("\r\n").split("\t")
    header2 = handle.readline().rstrip("\r\n").split("\t")
    header1 = [x.strip() for x in header1]
    header2 = [x.strip() for x in header2]
    assert len(header1) >= len(header2)
    header2 = header2 + [""]*(len(header1)-len(header2))
    header = ["%s %s" % (x1, x2) for (x1, x2) in zip(header1, header2)]
    header = [x.strip() for x in header]
    yield "\t".join(header)+"\n"

    x = handle.readline().strip()
    assert x == "-"*len(x)

    for cols in filelib.read_cols(handle):
        assert len(cols) == len(header), "%d %d" % (len(header), len(cols))
        yield "\t".join(cols)+"\n"

def parse_blat_psl(handle):
    # Yields hits as objects.  See blat for the members of the
    # objects.
    import filelib
    
    for d in filelib.read_row(_clean_blat_psl(handle), header=1):
        yield d

def calc_score_psl(psl):
    # Return perc_identity, score
    perc_identity = 100.0 - pslCalcMilliBad(psl, True)*0.1
    score = pslScore(psl)
    return perc_identity, score

# http://genome.ucsc.edu/FAQ/FAQblat
def pslCalcMilliBad(psl, isMrna):
    import math
    
    psl.match = int(psl.match)
    psl.rep__match = int(psl.rep__match)
    psl.mis__match = int(psl.mis__match)
    psl.Q_start, psl.Q_end = int(psl.Q_start), int(psl.Q_end)
    psl.T_start, psl.T_end = int(psl.T_start), int(psl.T_end)
    psl.Q_gap_count = int(psl.Q_gap_count)
    psl.T_gap_count = int(psl.T_gap_count)

    sizeMul = 1
    if pslIsProtein(psl):
        sizeMul = 3
    milliBad = 0

    qAliSize = sizeMul * (psl.Q_end - psl.Q_start)
    tAliSize = psl.T_end - psl.T_start
    aliSize = min(qAliSize, tAliSize)
    if aliSize <= 0:
        return 0
    sizeDif = qAliSize - tAliSize
    if sizeDif < 0:
        sizeDif = -sizeDif
        if isMrna:
            sizeDif = 0
    insertFactor = psl.Q_gap_count
    if not isMrna:
        insertFactor += psl.T_gap_count

    total = sizeMul * (psl.match+psl.rep__match+psl.mis__match)
    if total:
        milliBad = (psl.mis__match*sizeMul + insertFactor +
                    round(3*math.log(1+sizeDif)))
        milliBad = int(float(1000 * milliBad) / total)
    return milliBad

def pslIsProtein(psl):
    psl.block_count = int(psl.block_count)
    if type(psl.blockSizes) is type(""):
        x = [int(x) for x in psl.blockSizes.split(",") if x]
        psl.blockSizes = x
    if type(psl.qStarts) is type(""):
        x = [int(x) for x in psl.qStarts.split(",") if x]
        psl.qStarts = x
    if type(psl.tStarts) is type(""):
        x = [int(x) for x in psl.tStarts.split(",") if x]
        psl.tStarts = x

    psl.T_start, psl.T_end = int(psl.T_start), int(psl.T_end)
    psl.T_size = int(psl.T_size)

    lastBlock = psl.block_count - 1

    return (
        ((psl.strand == "+") and
         (psl.T_end == psl.tStarts[lastBlock] + 3*psl.blockSizes[lastBlock]))
        or
        ((psl.strand == "-") and
         (psl.T_start == (
        psl.T_size-(psl.tStarts[lastBlock] + 3*psl.blockSizes[lastBlock]))))
        )

def pslScore(psl):
    sizeMul = 1
    if pslIsProtein(psl):
        sizeMul = 3

    psl.match = int(psl.match)
    psl.rep__match = int(psl.rep__match)
    psl.mis__match = int(psl.mis__match)
    psl.Q_gap_count = int(psl.Q_gap_count)
    psl.T_gap_count = int(psl.T_gap_count)

    x = (sizeMul * (psl.match + (psl.rep__match/2)) -
         sizeMul * psl.mis__match - psl.Q_gap_count - psl.T_gap_count)
    return x
