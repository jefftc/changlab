from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_filename):
        from genomicode import filelib
        from genomicode import vcflib

        in_mvcf_node, back_mvcf_node = antecedents
        filelib.assert_exists_nz(in_mvcf_node.identifier)
        filelib.assert_exists_nz(back_mvcf_node.identifier)
        #print in_mvcf_node.identifier
        #print back_mvcf_node.identifier
        in_vcf = vcflib.read(in_mvcf_node.identifier)
        back_vcf = vcflib.read(back_mvcf_node.identifier)

        assert sorted(in_vcf.samples) == sorted(back_vcf.samples), \
               "Mismatched samples"

        # Make sure there are no duplicate sample names.
        x = {}.fromkeys(in_vcf.samples).keys()
        assert len(in_vcf.samples) == len(x), "Duplicate samples"

        # Parse out the read counts from the backfill vcf.
        sample2infolist = {}  # sample -> list of Info objects
        for sample in back_vcf.samples:
            x = [
                vcflib.parse_info(back_vcf, sample, i) 
                for i in range(back_vcf.matrix.num_annots())]
            sample2infolist[sample] = x

        all_chrom = back_vcf.matrix["#CHROM"]
        all_pos = [int(x) for x in back_vcf.matrix["POS"]]
        all_ref = back_vcf.matrix["REF"]
        all_alt = back_vcf.matrix["ALT"]
        sample2pos2info = {}  # sample -> (chrom, pos) -> Info
        for sample in sample2infolist:
            infolist = sample2infolist[sample]
            assert len(infolist) == len(all_chrom)
            assert len(infolist) == len(all_pos)
            pos2info = {}
            for x in zip(all_chrom, all_pos, all_ref, all_alt, infolist):
                chrom, pos, ref, alt, info = x
                # Skip this if no information.
                if info.num_ref is None and info.num_alt is None and \
                       info.total_reads is None and info.vaf is None:
                    continue
                pos2info[(chrom, pos, ref, alt)] = info
            sample2pos2info[sample] = pos2info

        # Find pos in final VCF file that matches the backfill VCF
        # file.
        all_chrom = in_vcf.matrix["#CHROM"]
        all_pos = [int(x) for x in in_vcf.matrix["POS"]]
        all_ref = in_vcf.matrix["REF"]
        all_alt = in_vcf.matrix["ALT"]
        matches = []  # List of (sample, variant_num, backfill_info).
        for sample in in_vcf.samples:
            pos2info = sample2pos2info[sample]
            for i, x in enumerate(zip(all_chrom, all_pos, all_ref, all_alt)):
                chrom, pos, ref, alt = x
                # Skip if there is no backfill information.
                bf_info = pos2info.get((chrom, pos, ref, alt))
                if bf_info is None:
                    continue
                # Skip if I already have information.
                x = vcflib.parse_info(in_vcf, sample, i)
                if not(x.num_ref is None and x.num_alt is None and 
                       x.total_reads is None and x.vaf is None):
                    continue
                
                x = sample, i, bf_info
                matches.append(x)

        # Update the read counts from annotated VCF file.
        for x in matches:
            sample, var_num, fill_info = x
            vcflib.update_info(in_vcf, sample, var_num, fill_info)

        vcflib.write(out_filename, in_vcf)


    def name_outfile(self, antecedents, user_options):
        return "backfilled.vcf"
