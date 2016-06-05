from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_filename):
        import shutil
        from genomicode import filelib
        from Betsy import module_utils as mlib

        simple_file = in_data.identifier
        metadata = {}

        x = mlib.get_user_option(user_options, "remove_sample")
        x = x.split(",")
        x = [x.strip() for x in x]
        remove_samples = x

        x = mlib.get_user_option(
            user_options, "apply_filter", allowed_values=["no", "yes"])
        apply_filter = (x == "yes")

        min_reads = mlib.get_user_option(
            user_options, "filter_by_min_reads", not_empty=True, type=int)
        assert min_reads >= 0 and min_reads < 10000
        
        CALLER2FILTERFN = {
            "MuTect" : is_pass_mutect,
            "VarScan2" : is_pass_varscan,
            "Strelka" : is_pass_strelka,
            "SomaticSniper" : is_pass_somaticsniper,
            "mutationSeq" : is_pass_jointsnvmix,
            }

        TEMPFILE = "temp.txt"
        handle = open(TEMPFILE, 'w')
        it = filelib.read_row(simple_file, header=1)
        print >>handle, "\t".join(it._header)
        for d in it:
            # remove_sample
            if d.Sample in remove_samples:
                continue

            # apply_filter
            if apply_filter:
                assert d.Caller in CALLER2FILTERFN, \
                       "Unknown caller: %s" % d.Caller
                is_pass_fn = CALLER2FILTERFN[d.Caller]
                if not is_pass_fn(d):
                    continue

            # filter_by_min_reads
            if min_reads > 0:
                x = d.Total_Reads
                x = x.split(",")
                x = [int(x) for x in x]
                x = max(x)
                total_reads = x
                if total_reads < min_reads:
                    continue
            print >>handle, "\t".join(d._cols)
        handle.close()

        shutil.move(TEMPFILE, out_filename)

        return metadata
            
    
    def name_outfile(self, antecedents, user_options):
        return "variants.txt"


def is_pass_mutect(d):
    # FILTER:
    # PASS
    # REJECT

    # Keep FILTER==PASS
    assert d.Filter in ["PASS", "REJECT"]
    if d.Filter == "REJECT":
        return False
    return True


def is_pass_varscan(d):
    # FILTER:
    # PASS
    # indelError

    # Keep FILTER==PASS
    assert d.Filter in ["PASS", "indelError"]
    if d.Filter == "indelError":
        return False
    return True


def is_pass_strelka(d):
    # FILTER:
    # PASS
    # BCNoise
    # BCNoise;DP
    # BCNoise;QSS_ref
    # BCNoise;QSS_ref;DP
    # BCNoise;QSS_ref;SpanDel
    # DP
    # DP;SpanDel
    # QSS_ref
    # QSS_ref;DP
    # QSS_ref;DP;SpanDel
    # QSS_ref;SpanDel

    # While VCF file is semicolon separated, SimpleVariantFile is
    # comma separated.

    assert d.Filter.find(";") < 0
    filters = d.Filter.split(",")
    if len(filters) > 1:
        return False
    if filters[0] != "PASS":
        return False
    return True


def is_pass_somaticsniper(d):
    # FILTER:
    # .
    return True


def is_pass_jointsnvmix(d):
    # FILTER:
    # PASS
    return True
