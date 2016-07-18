from Module import AbstractModule

# Columns:
# Caller
# Sample
# Chrom
# Pos
# Ref
# Alt
# Source         DNA or RNA
# Num Ref
# Num Alt
# Total Reads
# VAF
# Filter
# Call
# GQ

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_filename):
        import shutil
        from genomicode import parallel

        vcf_folder = in_data
        vcf_files = find_vcf_files(vcf_folder.identifier)
        metadata = {}

        TEMPFILE = "temp.txt"
        handle = open(TEMPFILE, 'w')
        header = (
            "Caller", "File", "Sample", "Chrom", "Pos", "Ref", "Alt", "Source",
            "Num Ref", "Num Alt", "Total Reads", "VAF", "Filter", "Call", "GQ")
        print >>handle, "\t".join(header)
        handle.close()

        # Write out data from each of the VCF files.
        jobs = []
        for x in vcf_files:
            filestem, filename = x

            # filestem   197B-MG
            # filename   /data/jchang/biocore/call01/radia.vcf/197B-MG.vcf

            args = filename, filestem, header, TEMPFILE
            x = summarize_vcf_file, args, {}
            jobs.append(x)
        parallel.pyfun(jobs, num_procs=num_cores, lock_keyword="lock")
        metadata["num_cores"] = num_cores

        shutil.move(TEMPFILE, out_filename)

        return metadata
            
    
    def name_outfile(self, antecedents, user_options):
        return "variants.txt"


def summarize_vcf_file(filename, filestem, header, outfilename, lock):
    from genomicode import hashlib
    from genomicode import vcflib

    vcf = vcflib.read(filename)

    lines = []
    for i in range(vcf.num_variants()):
        var = vcflib.get_variant(vcf, i)

        caller_name = var.caller.name
        ref = var.ref
        alt = ",".join(var.alt)
        filter_str = vcf.caller.get_filter(var)

        for sample in var.samples:
            # If sample begins with an integer, there may be a
            # "X" pre-pended to it.  Try to detect this case
            # and fix it.
            clean_sample = sample
            if sample == hashlib.hash_var(filestem):
                clean_sample = filestem

            source = "DNA"
            if caller_name == "Radia":
                # DNA    <clean_sample>       196B-lung
                # RNA    <clean_sample>_RNA   196B-lung_RNA
                # Figure out whether this is RNA and fix it.
                if clean_sample.endswith("_RNA"):
                    clean_sample = clean_sample[:-4]
                    source = "RNA"

            genodict = var.sample2genodict[sample]
            call = vcflib.get_call(var, sample)

            num_ref = vcflib._format_vcf_value(call.num_ref, None_char="")
            num_alt = vcflib._format_vcf_value(call.num_alt, None_char="")
            total_reads = vcflib._format_vcf_value(
                call.total_reads, None_char="")
            vaf = vcflib._format_vcf_value(call.vaf, None_char="")
            call_str = vcflib._format_vcf_value(call.call, None_char="")
            GQ = genodict.get("GQ", "")
            if GQ in [None, "."]:
                GQ = ""

            x = caller_name, filestem, clean_sample, var.chrom, var.pos, \
                ref, alt, source, \
                num_ref, num_alt, total_reads, vaf, filter_str, call_str, GQ
            assert len(x) == len(header)
            x = "\t".join(map(str, x))
            lines.append(x)

            if len(lines) >= 100000:
                x = "\n".join(lines)+"\n"
                lock.acquire()
                handle = open(outfilename, 'a')
                handle.write(x)
                handle.close()
                lock.release()
                lines = []

    x = "\n".join(lines)+"\n"
    lock.acquire()
    handle = open(outfilename, 'a')
    handle.write(x)
    handle.close()
    lock.release()
    

def find_vcf_files(vcf_path):
    # Return list of (<sample>, <filename>).
    
    import os
    from genomicode import filelib
    from genomicode import vcflib

    filenames = filelib.list_files_in_path(
        vcf_path, endswith=".vcf", case_insensitive=True)
    
    # Format:
    # <path>/<sample>.vcf
    vcf_files = []
    for filename in filenames:
        p, f = os.path.split(filename)
        sample = os.path.splitext(f)[0]
        #caller = vcflib.identify_caller(filename)
        #assert caller is not None, "Unknown caller: %s" % filename
        x = sample, filename
        vcf_files.append(x)
    return vcf_files
