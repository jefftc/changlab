from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_filename):
        import os
        import stat
        import shutil
        from genomicode import filelib
        from genomicode import parallel
        from Betsy import module_utils as mlib

        vcf_node = in_data
        # Some callers, like jointsnvmix, will create vcf files for
        # each chromosome.  To avoid picking these up, only accept
        # .vcf files from the top level.
        vcf_filenames = filelib.list_files_in_path(
            vcf_node.identifier, endswith=".vcf", toplevel_only=True)
        assert vcf_filenames, "No .vcf files."
        metadata = {}

        # TODO: find the version number of these tools.
        bgzip = mlib.findbin("bgzip")
        tabix = mlib.findbin("tabix")
        bcftools = mlib.findbin("bcftools")
        sq = parallel.quote
        
        tmp_path = "indexed.vcf"
        tmp_path = os.path.realpath(tmp_path)
        filelib.safe_mkdir(tmp_path)
        

        # list of (sample, in_filename, tmp_filename)
        jobs = []
        for in_filename in vcf_filenames:
            p, f = os.path.split(in_filename)
            sample = os.path.splitext(f)[0]
            tmp_filename = os.path.join(tmp_path, f)
            x = sample, in_filename, tmp_filename
            jobs.append(x)
        jobs.sort()   # Merge them in order of sample

        # Ignore VCF files that don't have any variants.
        jobs2 = []
        for x in jobs:
            sample, in_filename, tmp_filename = x
            if os.stat(in_filename)[stat.ST_SIZE] <= 0:
                continue
            jobs2.append(x)
        jobs = jobs2

        # If there are no VCF files with any variants, then just
        # create an empty outfile and return.
        if not jobs:
            open(out_filename, 'w')
            return
            

        # Copy all the VCF files to a temporary directory.
        for x in jobs:
            sample, in_filename, tmp_filename = x
            shutil.copy2(in_filename, tmp_filename)

        # NextGENe creates broken VCF files.  Fix them.
        for x in jobs:
            sample, in_filename, tmp_filename = x
            fix_nextgene_vcf(tmp_filename)

        ## # Since we are merging the files, we need to make sure that
        ## # each file has a unique name.  If the names aren't unique,
        ## # then make them unique by adding the name of the file.
        ## all_unique = True
        ## seen = {}
        ## for x in jobs:
        ##     sample, in_filename, tmp_filename = x
        ##     samples = _get_samples_from_vcf(tmp_filename)
        ##     for s in samples:
        ##         if s in seen:
        ##             all_unique = False
        ##             break
        ##         seen[s] = 1
        ##     if not all_unique:
        ##         break
        ## if not all_unique:
        ##     for x in jobs:
        ##         sample, in_filename, tmp_filename = x
        ##         _uniquify_samples_in_vcf(tmp_filename, sample)

        # Compress the VCF files.
        # bgzip file.vcf
        commands = []
        for x in jobs:
            sample, in_filename, tmp_filename = x
            x = "%s %s" % (sq(bgzip), sq(tmp_filename))
            commands.append(x)
        parallel.pshell(commands, max_procs=num_cores, path=tmp_path)
        x = ["%s.gz" % x[-1] for x in jobs]
        filelib.assert_exists_nz_many(x)

        # Index the VCF files.
        # tabix -p vcf file.vcf.gz
        commands = []
        for x in jobs:
            sample, in_filename, tmp_filename = x
            x = "%s -p vcf %s.gz" % (sq(tabix), sq(tmp_filename))
            commands.append(x)
        parallel.pshell(commands, max_procs=num_cores, path=tmp_path)
        x = ["%s.gz.tbi" % x[-1] for x in jobs]
        filelib.assert_exists_nz_many(x)
        
        # Run bcftools
        ## For VCF files from somatic calls, the germline sample will
        ## be duplicated.  Add --force-samples to make sure this is
        ## still merged.

        # Since we need to append all the VCF files, it's easy to run
        # into error:
        # OSError: [Errno 7] Argument list too long
        # 
        # To reduce the chance of this, figure out the path of the
        # tmp_filename, and run the analysis in that path so we can
        # use relative filenames.
        tmp_path = None
        for x in jobs:
            sample, in_filename, tmp_filename = x
            path, file_ = os.path.split(tmp_filename)
            if tmp_path is None:
                tmp_path = path
            assert path == tmp_path
        
        cmd = [
            sq(bcftools),
            "merge", 
            "-o %s" % sq(out_filename),
            "-O v",
            "--force-samples", 
            ]
        for x in jobs:
            sample, in_filename, tmp_filename = x
            path, file_ = os.path.split(tmp_filename)
            assert path == tmp_path
            cmd.append("%s.gz" % file_)
        x = " ".join(cmd)
        parallel.sshell(x, path=tmp_path)
        filelib.assert_exists_nz(out_filename)

        return metadata


    def name_outfile(self, antecedents, user_options):
        return "merged.vcf"


def _get_samples_from_vcf(filename):
    lines, header_i, samples = _read_vcf(filename)
    return samples


def _read_vcf(filename):
    # Return a tuple of:
    # - a list of lines.  Each line is a list of columns.
    # - the index of the header row (or None)
    # - the sample names
    from genomicode import filelib
    
    lines = [x for x in filelib.read_cols(filename)]
    header_i = None
    for i, cols in enumerate(lines):
        if cols[0] == "#CHROM":
            header_i = i
            break
    assert header_i is not None, "Could not find #CHROM: %s" % filename

    header = lines[header_i]
    x = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
         "INFO", "FORMAT"]
    assert header[:len(x)] == x, "Unknown format: %s" % header
    samples = header[len(x):]
    return lines, header_i, samples

    
def _uniquify_samples_in_vcf(filename, to_append):
    from genomicode import hashlib
    
    lines, header_i, samples = _read_vcf(filename)

    header = lines[header_i]
    header1 = header[:-len(samples)]
    header2 = header[-len(samples):]

    x = header2
    x = ["%s %s" % (x, to_append) for x in x]
    x = [hashlib.hash_var(x) for x in x]
    header2 = x
    lines[header_i] = header1 + header2

    handle = open(filename, 'w')
    for x in lines:
        print >>handle, "\t".join(x)
    handle.close()


def fix_nextgene_vcf(filename):
    # NextGENe adds a line that causes problems.  Get rid of it.
    # ##sgmutationstatistics=<SampleName=10_S16.pjt,TotalCount=1688,
    # Homozygous=1113,Heterozygous=575,Substitutions=1420,Insertions=129,
    # Deletions=148,Ti/Tv=1.39,Confirmed=0,Deleted=0,Added=0>
    lines = open(filename).readlines()
    x = [x for x in lines if not x.startswith("##sgmutationstatistics=")]
    if len(x) != len(lines):
        open(filename, 'w').writelines(x)
