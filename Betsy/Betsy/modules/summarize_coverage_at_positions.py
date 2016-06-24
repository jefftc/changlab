from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options,
        num_cores, outfile):
        from genomicode import filelib
        from genomicode import vcflib
        from Betsy import module_utils as mlib

        vcf_node = in_data
        vcf_filenames = filelib.list_files_in_path(
            vcf_node.identifier, endswith=".vcf", toplevel_only=True)
        assert vcf_filenames, "No .vcf files."
        metadata = {}

        jobs = []  # list of (filestem, vcf_filename)
        for vcf_filename in vcf_filenames:
            path, root, ext = mlib.splitpath(vcf_filename)
            x = root, vcf_filename
            jobs.append(x)

        vcf_objects = []
        for x in jobs:
            filestem, vcf_filename = x
            vcf = vcflib.read(vcf_filename)
            vcf_objects.append(vcf)

        # Make a list of all the samples.  Make sure they are unique.
        all_samples = []
        for vcf in vcf_objects:
            for sample in vcf.samples:
                assert sample not in all_samples
                all_samples.append(sample)

        # Make a list of all the positions.
        # sample name -> (chrom, pos) -> call
        sample2coord2call = {}
        for vcf in vcf_objects:
            for i in range(vcf.num_variants()):
                var = vcflib.get_variant(vcf, i)
                
                for sample in vcf.samples:
                    coord2call = sample2coord2call.get(sample, {})
                    coord = var.chrom, var.pos
                    call = vcflib.get_call(var, sample)
                    assert coord not in coord2call
                    coord2call[coord] = call
                    sample2coord2call[sample] = coord2call

        all_coord = {}
        for coord2call in sample2coord2call.itervalues():
            for coord in coord2call.iterkeys():
                all_coord[coord] = 1
        all_coord = sorted(all_coord)


        handle = open(outfile, 'w')
        header = ["Chrom", "Pos"] + all_samples
        print >>handle, "\t".join(header)
        for coord in all_coord:
            chrom, pos = coord
            pos_f = vcflib._format_vcf_value(pos)

            coverage = [""] * len(all_samples)
            for i, sample in enumerate(all_samples):
                coord2call = sample2coord2call.get(sample, {})
                call = coord2call.get(coord, None)
                if not call:
                    continue
                coverage[i] = vcflib._format_vcf_value(
                    call.total_reads, char_for_None="")
            
            x = [chrom, pos_f] + coverage
            assert len(x) == len(header)
            print >>handle, "\t".join(map(str, x))
        handle.close()
            
        return metadata

        
    def name_outfile(self, antecedents, user_options):
        return "coverage.txt"
