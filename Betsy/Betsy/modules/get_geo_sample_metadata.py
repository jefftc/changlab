"""

Methods:
run
name_outfile
make_unique_hash
get_out_attributes
find_antecedents

"""


def run(network, antecedents, out_attributes, user_options, num_cores):
    in_data = antecedents

    from genomicode import geolib
    from genomicode import jmath

    from Betsy import module_utils
    from Betsy import bie3
    from Betsy import rulebase

    convert_NA = user_options.get("set_NA_to")

    # Input should be a GEOSeriesMatrixFile.
    idata_in = antecedents
    filename = idata_in.identifier
    outfile = name_outfile(in_data, user_options)

    # Get the sample data and write it out.
    matrix = geolib._extract_sm_sample_meta(filename)
    matrix = geolib._clean_sm_sample_meta(matrix)
    matrix = geolib._prettify_sm_sample_meta(matrix)
    matrix = jmath.transpose(matrix)   # each column is an annotation
    if convert_NA != "NA":
        for i in range(1, len(matrix)):
            for j in range(len(matrix[i])):
                # Do case sensitive?
                if matrix[i][j] == "NA":
                    matrix[i][j] = convert_NA
    outhandle = open(outfile, 'w')
    for x in matrix:
        print >>outhandle, "\t".join(map(str, x))
    outhandle.close()

    out_node = bie3.Data(rulebase.GEOSampleMetadata, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


# Not sure exactly what this does.  Looks like it makes a unique name
# for the output of this module.
def name_outfile(antecedents, user_options):
    import os

    # BUG: what if there are multiple GPLIDs?
    #original_file = module_utils.get_inputid(user_options['GSEID'])
    #filename = "download_geo_seriesmatrix_%s" % original_file
    filename = "sample_metadata.tsv"
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    from Betsy import module_utils
    
    identifier = antecedents.identifier
    return module_utils.make_unique_hash(
        identifier, pipeline, out_attributes, user_options)


def set_out_attributes(antecedents, out_attributes):
    return out_attributes


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    from Betsy import module_utils
    
    data_node = module_utils.find_antecedents(
        network, module_id, user_attributes, pool)
    return data_node
