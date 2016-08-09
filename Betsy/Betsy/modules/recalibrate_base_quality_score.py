from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import filelib
        from genomicode import parallel
        from genomicode import alignlib
        from Betsy import module_utils

        bam_node, ref_node, report_node = antecedents
        bam_filenames = module_utils.find_bam_files(bam_node.identifier)
        assert bam_filenames, "No .bam files."
        report_filenames = filelib.list_files_in_path(
            report_node.identifier, endswith=".grp")
        assert report_filenames, "No .grp files."
        ref = alignlib.create_reference_genome(ref_node.identifier)
        filelib.safe_mkdir(out_path)
        metadata = {}

        assert len(bam_filenames) == len(report_filenames), \
               "Should have a .grp file for each .bam file."
        sample2bamfilename = {}
        for filename in bam_filenames:
            p, f = os.path.split(filename)
            sample, ext = os.path.splitext(f)
            assert sample not in sample2bamfilename
            sample2bamfilename[sample] = filename
        sample2reportfilename = {}
        for filename in report_filenames:
            p, f = os.path.split(filename)
            sample, ext = os.path.splitext(f)
            assert sample not in sample2reportfilename
            sample2reportfilename[sample] = filename
        assert len(sample2bamfilename) == len(sample2reportfilename)

        missing = [
            x for x in sample2bamfilename if x not in sample2reportfilename]
        assert not missing, "Missing grp files for %d bam files." % \
               len(missing)

        jobs = []  # list of (bam_filename, report_filename, out_filename)
        for sample in sample2bamfilename:
            bam_filename = sample2bamfilename[sample]
            report_filename = sample2reportfilename[sample]
            
            p, f = os.path.split(bam_filename)
            sample, ext = os.path.splitext(f)
            out_filename = os.path.join(out_path, "%s.bam" % sample)
            x = bam_filename, report_filename, out_filename
            jobs.append(x)

        # java -jar GenomeAnalysisTK.jar \
        #   -T PrintReads \
        #   -R reference.fasta \
        #   -BQSR recalibration_report.grp \
        #   -I input.bam -o output.bam

        # Make a list of commands.
        commands = []
        for x in jobs:
            in_filename, report_filename, out_filename = x
            nc = max(1, num_cores/len(jobs))
            x = alignlib.make_GATK_command(
                T="PrintReads", R=ref.fasta_file_full, nct=str(nc),
                BQSR=report_filename, I=in_filename, o=out_filename)
            commands.append(x)

        parallel.pshell(commands, max_procs=num_cores)
        metadata["num_cores"] = num_cores
        metadata["commands"] = commands

        # Make sure the analysis completed successfully.
        out_filenames = [x[-1] for x in jobs]
        filelib.assert_exists_nz_many(out_filenames)

        return metadata


    def name_outfile(self, antecedents, user_options):
        return "recalibrated.bam"

