from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        import shutil
        from genomicode import parallel
        from genomicode import filelib
        from genomicode import alignlib
        from Betsy import module_utils as mlib

        bam_filenames = mlib.find_bam_files(in_data.identifier)
        filelib.safe_mkdir(out_path)

        metadata = {}
        metadata["tool"] = "bam2fastx (unknown version)"

        # Somehow bam2fastx doesn't work if there are spaces in the
        # filename.  Make a temporary filename with no spaces, and
        # then rename it later.
        # Actually, may not be bam2fastx's fault.

        jobs = []
        for i, bam_filename in enumerate(bam_filenames):
            p, f, e = mlib.splitpath(bam_filename)
            #bai_filename = alignlib.find_bai_file(bam_filename)
            #assert bai_filename, "Missing index for: %s" % bam_filename
            #temp_bam_filename = "%d.bam" % i
            #temp_bai_filename = "%d.bam.bai" % i
            #temp_fa_filename = "%d.fa" % i
            fa_filename = os.path.join(out_path, "%s.fa" % f)
            x = filelib.GenericObject(
                bam_filename=bam_filename,
                #bai_filename=bai_filename,
                #temp_bam_filename=temp_bam_filename,
                #temp_bai_filename=temp_bai_filename,
                #temp_fa_filename=temp_fa_filename,
                fa_filename=fa_filename)
            jobs.append(x)
        bam2fastx = mlib.findbin("bam2fastx")

        # Link all the bam files.
        #for j in jobs:
        #    assert not os.path.exists(j.temp_bam_filename)
        #    #assert not os.path.exists(j.temp_bai_filename)
        #    os.symlink(j.bam_filename, j.temp_bam_filename)
        #    #os.symlink(j.bai_filename, j.temp_bai_filename)

        commands = []
        for j in jobs:
            # bam2fastx -A --fasta -o rqc14.fa rqc11.bam
            x = [
                mlib.sq(bam2fastx),
                "-A",
                "--fasta",
                #"-o", mlib.sq(j.temp_fa_filename),
                #mlib.sq(j.temp_bam_filename),
                "-o", mlib.sq(j.fa_filename),
                mlib.sq(j.bam_filename),
                ]
            x = " ".join(x)
            commands.append(x)
        metadata["commands"] = commands
        metadata["num_cores"] = num_cores
        parallel.pshell(commands, max_procs=num_cores)

        #for j in jobs:
        #    # Move the temporary files to the final location.
        #    shutil.move(j.temp_fa_filename, j.fa_filename)
        #    # Remove the link to the BAM file.
        #    os.unlink(j.temp_bam_filename)
        
        x = [j.fa_filename for x in jobs]
        filelib.assert_exists_nz_many(x)

        return metadata


    def name_outfile(self, antecedents, user_options):
        return "fasta"
