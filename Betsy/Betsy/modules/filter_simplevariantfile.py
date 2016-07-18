from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_filename):
        import shutil
        from genomicode import filelib
        from genomicode import vcflib
        from Betsy import module_utils as mlib

        simple_file = in_data.identifier
        metadata = {}

        x = mlib.get_user_option(user_options, "remove_samples")
        x = x.split(",")
        x = [x.strip() for x in x]
        remove_samples = x

        #x = mlib.get_user_option(
        #    user_options, "remove_radia_rna_samples",
        #    allowed_values=["no", "yes"])
        #remove_radia_rna_samples = (x == "yes")


        x = mlib.get_user_option(
            user_options, "apply_filter", allowed_values=["no", "yes"])
        apply_filter = (x == "yes")

        wgs_or_wes = mlib.get_user_option(
            user_options, "wgs_or_wes", not_empty=True,
            allowed_values=["wgs", "wes"])

        TEMPFILE = "temp.txt"
        handle = open(TEMPFILE, 'w')
        it = filelib.read_row(simple_file, header=1)
        print >>handle, "\t".join(it._header)
        for d in it:

            # Find the caller.
            for caller in vcflib.CALLERS:
                caller = caller()
                if caller.name == d.Caller:
                    break
            else:
                raise AssertionError, "Unknown caller: %s" % caller.name
            
            # remove_sample
            if d.Sample in remove_samples:
                continue
            #if remove_radia_rna_samples and d.Sample.endswith("_RNA"):
            #    continue

            # apply_filter
            if apply_filter:
                args = d.Filter,
                if d.Caller == "MuSE":
                    args = d.Filter, wgs_or_wes
                if not caller.is_pass(*args):
                    continue

            ## # filter_by_min_alt_reads
            ## if min_alt_reads > 0:
            ##     alt_reads = 0
            ##     if d.Num_Alt:
            ##         x = d.Num_Alt
            ##         x = x.split(",")
            ##         x = [int(x) for x in x]
            ##         x = max(x)
            ##         alt_reads = x
            ##     if alt_reads < min_alt_reads:
            ##         continue

            ## # filter_by_min_total_reads
            ## if min_total_reads > 0:
            ##     x = d.Total_Reads
            ##     x = x.split(",")
            ##     x = [int(x) for x in x]
            ##     x = max(x)
            ##     total_reads = x
            ##     if total_reads < min_total_reads:
            ##         continue
                
            ## # filter_by_min_GQ
            ## # If GQ isn't provided, then ignore this filter.
            ## if min_gq > 0 and d.GQ:
            ##     GQ = float(d.GQ)
            ##     if GQ < min_gq:
            ##         continue
            print >>handle, "\t".join(d._cols)
        handle.close()

        shutil.move(TEMPFILE, out_filename)

        return metadata
            
    
    def name_outfile(self, antecedents, user_options):
        return "variants.txt"
