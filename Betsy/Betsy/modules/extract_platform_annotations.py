from Module import AbstractModule


class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        from genomicode import filelib

        outhandle = open(outfile, 'w')
        extract_sf_platform(in_data.identifier, outhandle)
        outhandle.close()
        filelib.assert_exists_nz(outfile)

    def name_outfile(self, antecedents, user_options):
        return "platform.txt"



def extract_sf_platform(filename, outhandle):
    from genomicode import filelib
    
    handle = filelib.openfh(filename)
    while 1:
        line = handle.readline()
        if not line:
            raise AssertionError, "I could not find platform"
        #if line.startswith("^PLATFORM") and line.find(platform) >= 0:
        #    break
        # Assuming only one platform per file.
        if line.startswith("^PLATFORM"):
            break

    in_platform_table = 0
    for line in handle:
        if line.startswith("!platform_table_begin"):
            in_platform_table = 1
        elif line.startswith("!platform_table_end"):
            break
        elif in_platform_table:
            print >>outhandle, line,
    handle.close()
