#!/usr/local/bin/python

import sys, os

# read_signatures
# format_module
# form_get

SIGNATURE_PATH = "/opt/SIGNATURE"
sys.path = [SIGNATURE_PATH] + sys.path


def read_signatures(sigdb_path, desired_normalization, desired_tags):
    from genomicode import filefns
    
    opr = os.path.realpath
    opj = os.path.join

    filename = opj(sigdb_path, "signatures.txt")
    assert os.path.exists(filename), "Missing signatures.txt file."

    x = [x.upper() for x in desired_tags]
    desired_tags = {}.fromkeys(x)

    x = [x.upper() for x in desired_normalization]
    desired_normalization = {}.fromkeys(x)

    ds = []
    for d in filefns.read_row(filename, header=1):
        # Skip if not the right normalization.
        if d.Normalization.upper() not in desired_normalization:
            continue
        
        # Skip if the tags don't match.
        tags = [x.upper() for x in d.Tags.split()]
        for tag in tags:
            if tag in desired_tags:
                break
        else:
            # None of these tags matched any of the desired ones.
            continue
        
        # Skip if not all parameters supplied.
        if not d.Normalization:
            continue
        if not d.Genes or not d.Metagenes:
            continue
        if not d.Metagenes or not d.Quantile:
            continue
        if not d.Train0 or not d.Train1:
            continue

        # Find the training files.  If not found, then skip.
        train0 = opr(opj(sigdb_path, d.Train0))
        train1 = opr(opj(sigdb_path, d.Train1))
        if not os.path.exists(train0):
            train0 = train0 + ".gz"
        if not os.path.exists(train1):
            train1 = train1 + ".gz"
        if not os.path.exists(train0) or not os.path.exists(train1):
            continue
        d.Train0 = train0
        d.Train1 = train1
        
        # xls2txt converts all values to floats.
        d.xID = int(float(d.xID))
        d.Genes = int(float(d.Genes))
        d.Metagenes = int(float(d.Metagenes))
        ds.append(d)
    return ds

def format_module(form):
    from genomicode import genepatternfns as gp
    
    params = []

    # rma_expression_file
    default = form_get(form, "rma_expression_file_url", "")
    if not default:
        default = form_get(form, "rma_expression_file", "")
    x = gp.EFileChooser(
        "rma_expression_file", "RMA normalized expression data set. PCL, CDT, "
        "GCT, RES, ODF format.", [], default=default, optional=True)
    params.append(x)

    # mas5_expression_file
    default = form_get(form, "mas5_expression_file_url", "")
    if not default:
        default = form_get(form, "mas5_expression_file", "")
    x = gp.EFileChooser(
        "mas5_expression_file", "MAS5 normalized expression data set. PCL, "
        "CDT, GCT, RES, ODF format.", [], default=default, optional=True)
    params.append(x)

    # which_signatures: option to select signatures.
    which_signatures = form_get(form, "which_signatures", "all of them")
    x = gp.EDropdown(
        "which_signatures", "Which signatures do you want to run?",
        ["all of them", "I choose myself"], default=which_signatures,
        optional=True)
    params.append(x)

    # If the user wants to select_signatures themselves, then provide
    # an option for each signature.
    if which_signatures == "I choose myself":
        SIGDB_PATH = os.path.join(SIGNATURE_PATH, "sigdb")
        SIGNATURES = read_signatures(
            SIGDB_PATH, ["rma", "mas5"], ["production"])
        for sig in SIGNATURES:
            key = "sig_%s" % sig.Name
            default = form_get(form, key, "yes (default parameters)")
            description = (
                "Analyze the %s signature?  (Requires %s normalization.)" %
                (sig.Name, sig.Normalization))
            choices = [
                "yes (default parameters)", "yes (custom parameters)", "no"]
            x = gp.EDropdown(
                key, description, choices, default=default, optional=True)
            params.append(x)

            # If the user wants custom parameters for this signature,
            # then show the options for it.
            if default != "yes (custom parameters)":
                continue
            key_num_genes = "%s_num_genes" % key
            num_genes = form_get(form, key_num_genes, str(sig.Genes))
            x = gp.ETextBox(
                key_num_genes, "Number of genes.", default=num_genes,
                optional=True)
            params.append(x)

            key_num_metagenes = "%s_num_metagenes" % key
            num_metagenes = form_get(
                form, key_num_metagenes, str(sig.Metagenes))
            x = gp.ETextBox(
                key_num_metagenes, "Number of metagenes.",
                default=num_metagenes, optional=True)
            params.append(x)
            
            key_quantile = "%s_apply_quantile_normalization" % key
            default = "no"
            if sig.Quantile.upper() == "YES":
                default = "yes"
            quantile = form_get(form, key_quantile, default)
            x = gp.EDropdown(
                key_quantile, "Apply quantile normalization?", ["yes", "no"],
                default=quantile, optional=True)
            params.append(x)
            
            key_shiftscale = "%s_apply_shiftscale_normalization" % key
            default = "no"
            if sig.Shift_Scale.upper() == "YES":
                default = "yes"
            shiftscale = form_get(form, key_shiftscale, default)
            x = gp.EDropdown(
                key_shiftscale, "Apply shift-scale normalization?",
                ["yes", "no"], default=shiftscale, optional=True)
            params.append(x)

    params = [gp.elem2rtbparam(x) for x in params]
    x = gp.format_rtbparameter(params)
    return x

gp_imod_all_vars_dict = None
def _get_gp_imod_all_vars_dict(form):
    global gp_imod_all_vars_dict
    import urllib
    
    if gp_imod_all_vars_dict:
        return gp_imod_all_vars_dict
    GP_IMOD_ALL_VARS = "gp_imod_all_vars"
    if GP_IMOD_ALL_VARS not in form:
        return {}
    
    all_vars = {}
    value = form.getlist(GP_IMOD_ALL_VARS)[0]
    for x in value.split("&"):
        k, v = x.split("=", 1)
        # Unquote after we split on "&" and "=".
        k = urllib.unquote_plus(k)
        v = urllib.unquote_plus(v)
        all_vars[k] = v
    gp_imod_all_vars_dict = all_vars
    return gp_imod_all_vars_dict

def form_get(form, key, default):
    # If gp_imod_all_vars is given, it overrides anything else that
    # was specified.
    all_vars = _get_gp_imod_all_vars_dict(form)
    if key in all_vars:
        return all_vars[key]
    value = default
    if key in form:
        value = form.getlist(key)[0]
    return value

def main():
    import cgi
    
    print "Content-type: text/plain\n"

    # gp_modulename    name of module
    # <parameters>
    # ...
    form = cgi.FieldStorage()

    # Format the module as a list of RTBParameters.
    x = format_module(form)
    callback = form_get(form, "gp_callback", "")
    if callback:
        # JSONP request.
        x = "%s(%s)" % (callback, x)
    # Print out the JSON representation of the list of RTBParameters.
    print x,

    # Debug: write the parameters out to a debug file.
    if 0:
        # Need to make this file world-writeable so web server can
        # write to it.
        handle = open("/tmp/test.txt", 'w')
        print >>handle, "Parameters"
        for name in form:
            values = form.getlist(name)
            for value in values:
                print >>handle, "%s=%s" % (name, value)
        print >>handle
        print >>handle, "JSON parameters"
        print >>handle, x
        handle.close()

if __name__ == '__main__':
    main()
