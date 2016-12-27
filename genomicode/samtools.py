"""Code to work with samtools.

Functions:
parse_sam
reconstruct_sequence

Classes:
SamAlignment

"""
class SAMAlignment:
    def __init__(self, qname, flag, rname, pos, mapq, cigar, rnext, pnext,
                 tlen, seq, qual, tags):
        # tags should be dict of tag -> (type, value)
        assert type(flag) is type(0)
        assert type(pos) is type(0)
        assert type(mapq) is type(0)
        assert type(tlen) is type(0)
        assert type(tags) is type({})
        self.qname = qname
        self.flag = flag
        self.rname = rname
        self.pos = pos
        self.mapq = mapq
        self.cigar = cigar
        self.rnext = rnext
        self.pnext = pnext
        self.tlen = tlen
        self.seq = seq
        self.qual = qual
        self.tags = tags.copy()


def parse_sam(file_or_handle):
    # yield SAMAlignment objects
    from genomicode import filelib


    # Somehow, csv raises errors on some BAM files read directly with
    # "samtools view" (via the subprocess module).  Just implement our
    # own column splitting.
    #for cols in filelib.read_cols(file_or_handle):
    handle = filelib.openfh(file_or_handle)
    for line in handle:
        cols = line.rstrip("\r\n").split("\t")
        assert len(cols) >= 11, "Invalid line (%d):\n%s" % (len(cols), line)
        qname = cols[0]
        flag = int(cols[1])
        rname = cols[2]
        pos = int(cols[3])
        mapq = int(cols[4])
        cigar = cols[5]
        rnext = cols[6]
        pnext = cols[7]
        tlen = int(cols[8])
        seq = cols[9]
        qual = cols[10]
        tags = {}
        for i in range(11, len(cols)):
            x = cols[i]
            x = x.split(":", 2)
            assert len(x) == 3, cols[i]
            tag, type_, value = x
            assert tag not in tags
            tags[tag] = (type_, value)
        x = SAMAlignment(
            qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq,
            qual, tags)
        yield x


def _split_cigar(cigar_string):
    # Return a list of (number, operation).  If cigar_string is "*",
    # then return an empty list.
    import re

    if cigar_string == "*":
        return []

    parsed = []
    while cigar_string:
        m = re.match(r"\d+\w", cigar_string)
        assert m is not None
        x = m.group(0)
        cigar_string = cigar_string[len(x):]
        num, op = int(x[:-1]), x[-1]
        x = num, op
        parsed.append(x)
    return parsed


def reconstruct_sequence(seq, pos, cigar):
    # Figure out the sequence from a CIGAR string.
    # Return a dictionary of position -> base
    pos2base = {}
    for (n, op) in _split_cigar(cigar):
        if op == "S":       # soft clipping
            seq = seq[n:]
            #pos += n       # clip sequence, don't change position.
        elif op == "H":     # hard clipping
            pass            # not entirely sure we don't increment pos
        elif op == "I":
            assert n < 100, n
            assert n <= len(seq)
            p = (pos - 1) + 0.01
            for i in range(n):
                #if p not in pos2bases:
                #    pos2bases[p] = []
                #pos2bases[p].append(seq[0])
                assert p not in pos2base
                pos2base[p] = seq[0]
                seq = seq[1:]
                p += 0.01
        elif op == "M":
            assert n <= len(seq)
            for i in range(n):
                #if pos not in pos2bases:
                #    pos2bases[pos] = []
                #pos2bases[pos].append(seq[0])
                assert pos not in pos2base
                pos2base[pos] = seq[0]
                seq = seq[1:]
                pos += 1
        elif op == "D":
            pos += n
        else:
            raise AssertionError, "Unknown operation %s" % op
    return pos2base
