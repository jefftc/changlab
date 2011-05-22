"""

Functions:
primer3
primer3_core

parse

"""
import os, sys

def primer3(sequence, **params):
    # See primer3_core for more options.
    # Return list of (left_primer, right_primer, product_size)
    from StringIO import StringIO
    
    handle = StringIO()
    primer3_core(sequence, outhandle=handle, **params)
    handle.seek(0)
    #print handle.read(); sys.exit(0)
    return parse(handle)

def primer3_core(
    sequence, product_size=None, num_return=None, target=None, 
    format_output=None, primer3_core_app=None, outhandle=None, **params):
    # product_size  (75, 100) or [(100,125),(150,300)]  default 100-300
    # num_return    5           default 5
    # target        (25, 50)    target base 25 (1-based), length 50
    import subprocess
    import tempfile
    import config

    primer3_core_app = primer3_core_app or config.primer3_PRIMER3_BIN
    outhandle = outhandle or sys.stdout

    # Set the parameters to the webpage defaults.
    # 
    # Can reproduce results on web site if we set on the website:
    #   Targets
    #   Product Size Ranges
    # 
    # Default Product Size Range for web is:
    # 150-250 100-300 301-400 401-500 501-600 601-700 701-850 851-1000
    defaults = {
        "PRIMER_MAX_END_STABILITY" : 9.0,
        "PRIMER_MAX_TEMPLATE_MISPRIMING" : 12.00,
        "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING" : 24.00,
        "PRIMER_LIBERAL_BASE" : 1,
        "PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS" : 0,
        "PRIMER_QUALITY_RANGE_MAX" : 100,
        #"PRIMER_PAIR_WT_IO_PENALTY" : 1.0,
        }
    for name, value in defaults.iteritems():
        params[name] = params.get(name, value)

    cmd = [primer3_core_app]
    if format_output:
        cmd.append("-format_output")
    cmd.append("-strict_tags")
    cmd = " ".join(cmd)

    p = subprocess.Popen(
        cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    w, r = p.stdin, p.stdout
    w.write("SEQUENCE_ID=test\n")
    #w.write("SEQUENCE=%s\n" % sequence)  # obsolete
    w.write("SEQUENCE_TEMPLATE=%s\n" % sequence)
    if target is not None:
        assert len(target) == 2
        base, length = target
        assert base-1+length <= len(sequence)
        w.write("SEQUENCE_TARGET=%d,%d\n" % (base, length))
    if product_size:
        # Figure out if it's a tuple or a list of tuples.
        if type(product_size[0]) is type(0):
            assert len(product_size) == 2
            product_size = [product_size]
        for x in product_size:
            assert len(x) == 2
        # Format the product size.
        # space separated listed of <x>-<y>
        sizes = ["%d-%d" % x for x in product_size]
        x = " ".join(sizes)
        w.write("PRIMER_PRODUCT_SIZE_RANGE=%s\n" % x)
    if num_return is not None:
        assert num_return >= 1
        w.write("PRIMER_NUM_RETURN=%d\n" % num_return)
    for name, value in params.iteritems():
        w.write("%s=%s\n" % (name, value))

    w.write("=\n")
    w.close()

    for line in r:
        outhandle.write(line)

class PrimerInfo:
    # o pos is 0-based.
    # o The right primer is 5' to 3' for the primer.  Need to revcomp
    #   to compare against the original sequence.
    # o The position of right primer is relative to the 5' end of the
    #   primer (3' end of the sequence).  To get the position of the
    #   primer on a 0-based sequence, do:
    #   pos - length + 1
    def __init__(self, seq, pos, length, tm, gc_percent):
        self.seq = seq
        self.pos = pos
        self.length = length
        self.tm = tm
        self.gc_percent = gc_percent

def _parse_h(handle):
    # Yield: None or <num>, key, value
    import re
    
    for line in handle:
        assert not line.startswith("Unrecognized"), line.rstrip()
        assert "=" in line
        key, value = line.split("=", 1)
        key, value = key.strip(), value.strip()
        if not key and not value:
            break
        
        primer_id = None
        m = re.search(r"_(\d+)", key)
        if m:
            primer_id = int(m.group(1))
        yield primer_id, key, value

def parse(handle):
    # Return list of (left (PrimerInfo), right (PrimerInfo), product_size).
    import re

    # Which input parameters to ignore.
    #input_keys = [
    #    "SEQUENCE", "PRIMER_PRODUCT_SIZE_RANGE", "PRIMER_NUM_RETURN"]

    def data2obj(data, whichone):
        #print "START"
        #for k in data.keys():
        #    print k
        #print "END"
        seq = data["PRIMER_%s_SEQUENCE" % whichone]
        x = data["PRIMER_%s" % whichone]
        x, y = x.split(",")
        pos = int(x)
        length = int(y)
        tm = float(data["PRIMER_%s_TM" % whichone])
        gc_percent = float(data["PRIMER_%s_GC_PERCENT" % whichone])
        return PrimerInfo(seq, pos, length, tm, gc_percent)

    primers = []
    
    PRODUCT_SIZE = "PRIMER_PAIR_PRODUCT_SIZE"
    data = {}   # key -> value
    prev_primer_id = None
    for primer_id, key, value in _parse_h(handle):
        if primer_id != prev_primer_id:
            if prev_primer_id is not None:
                left_primer = data2obj(data, "LEFT")
                right_primer = data2obj(data, "RIGHT")
                product_size = int(data.get(PRODUCT_SIZE, 0))
                x = left_primer, right_primer, product_size
                primers.append(x)
            data = {}
        prev_primer_id = primer_id
        
        # PRIMER_PAIR_1_PENALTY ->  PRIMER_PAIR_PENALTY
        # PRIMER_LEFT_1         ->  PRIMER_LEFT
        key = re.sub(r"_\d+", "", key)
        data[key] = value
    if primer_id is not None:
        left_primer = data2obj(data, "LEFT")
        right_primer = data2obj(data, "RIGHT")
        product_size = int(data.get(PRODUCT_SIZE, 0))
        x = left_primer, right_primer, product_size
        primers.append(x)
    return primers
