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

        annot_header = ["Chrom", "Pos", "Ref", "Alt"]
        matrix = SimpleVariantMatrix.make_matrix(
            samples, callers, annot_header, coord_data, call_data)
        SimpleVariantMatrix.write(out_filename, matrix)

        return metadata

    
    def name_outfile(self, antecedents, user_options):
        return "calls.txt"
