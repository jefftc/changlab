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

        x = mlib.get_user_option(
            user_options, "apply_filter", allowed_values=["no", "yes"])
        apply_filter = (x == "yes")

        wgs_or_wes = mlib.get_user_option(
            user_options, "wgs_or_wes", not_empty=True,
            allowed_values=["wgs", "wes"])

        name2caller = {}  # name -> Caller object
        for caller in vcflib.CALLERS:
            caller = caller()
            assert caller.name not in name2caller
            name2caller[caller.name] = caller

        TEMPFILE = "temp.txt"
        handle = open(TEMPFILE, 'w')
        it = filelib.read_row(simple_file, header=1)
        print >>handle, "\t".join(it._header)
        for d in it:
            # Find the caller.
            assert d.Caller in name2caller, "Unknown caller: %s" % d.Caller
            caller = name2caller[d.Caller]
            
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

            print >>handle, "\t".join(d._cols)
        handle.close()

        shutil.move(TEMPFILE, out_filename)

        return metadata
            
    
    def name_outfile(self, antecedents, user_options):
        return "variants.txt"
