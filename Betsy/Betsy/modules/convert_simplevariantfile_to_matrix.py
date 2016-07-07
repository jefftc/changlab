from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_filename):
        from genomicode import filelib
        from genomicode import SimpleVariantMatrix

        simple_file = in_data.identifier
        metadata = {}

        # Read all in memory.  Hopefully, not too big.
        ds = []
        for d in filelib.read_row(simple_file, header=-1): 
            ds.append(d)
            #if len(ds) > 50000:  # DEBUG
            #    break

        # MuSE sometimes has alternates.
        # Alt       A,C
        # Num_Alt  13,0
        # VAF      0.19,0.0
        # Detect this and fix it.  Take the alternate with the highest VAF.
        for d in ds:
            if d.Num_Alt.find(",") < 0:
                continue
            x1 = d.Num_Alt.split(",")
            x2 = d.VAF.split(",")
            assert len(x1) == len(x2)
            x1 = map(int, x1)
            x2 = map(float, x2)
            max_vaf = max_i = None
            for i in range(len(x2)):
                if max_vaf is None or x2[i] > max_vaf:
                    max_vaf = x2[i]
                    max_i = i
            assert max_i is not None
            d.Num_Alt = str(x1[max_i])
            d.VAF = str(x2[max_i])

        # Make a list of all the positions.
        positions = {}  # (Chrom, Pos) -> 1
        for d in ds:
            positions[(d.Chrom, int(d.Pos))] = 1
        positions = sorted(positions)

        # Make a list of all the callers.
        callers = {}
        for d in ds:
            callers[d.Caller] = 1
        callers = sorted(callers)

        # Make a list of all the samples.
        samples = {}
        for d in ds:
            samples[d.Sample] = 1
        samples = sorted(samples)

        # Make a list of the coordinates.
        coord_data = {}
        for d in ds:
            x = d.Chrom, int(d.Pos), d.Ref, d.Alt
            coord_data[x] = 1
        coord_data = sorted(coord_data)

        call_data = []
        for d in ds:
            num_ref = num_alt = vaf = None
            if d.Num_Ref:
                num_ref = int(d.Num_Ref)
            if d.Num_Alt:
                num_alt = int(d.Num_Alt)
            if d.VAF:
                vaf = float(d.VAF)
            if num_ref is None and num_alt is None and vaf is None:
                continue
            call = SimpleVariantMatrix.Call(num_ref, num_alt, vaf)
            x = d.Chrom, int(d.Pos), d.Ref, d.Alt, d.Sample, d.Caller, call
            call_data.append(x)

        # Make a named matrix that shows the number of callers that
        # called a variant at each position for each sample.
        named_data = []
        samp2coord2nc = {}  # sample -> chrom, pos, ref, alt -> num_callers
        for x in call_data:
            chrom, pos, ref, alt, sample, caller, call = x
            coord = chrom, pos, ref, alt
            if sample not in samp2coord2nc:
                samp2coord2nc[sample] = {}
            nc = samp2coord2nc[sample].get(coord, 0) + 1
            samp2coord2nc[sample][coord] = nc
        headers = samples
        all_annots = []
        for sample in samples:
            annots = []
            for coord in coord_data:
                nc = samp2coord2nc.get(sample, {}).get(coord, "")
                annots.append(nc)
            all_annots.append(annots)
        if headers:
            x = "Num Callers", headers, all_annots
            named_data.append(x)

        annot_header = ["Chrom", "Pos", "Ref", "Alt"]
        matrix = SimpleVariantMatrix.make_matrix(
            samples, callers, annot_header, coord_data, named_data,
            call_data)
        SimpleVariantMatrix.write(out_filename, matrix)

        return metadata

    
    def name_outfile(self, antecedents, user_options):
        return "calls.txt"
