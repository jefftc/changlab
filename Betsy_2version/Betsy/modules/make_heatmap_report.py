#make_cluster_report.py

import os
import shutil
import imghdr
import time
from Betsy import bie
from Betsy import rulebase
from Betsy import config
from Betsy import module_utils
from Betsy import hash_method

def run(data_node,parameters, network):
    outfile=name_outfile(data_node)
    result_files = []
    filename = data_node.attributes['filename']
    new_name = os.path.split(filename)[-1]
    if os.path.isdir(filename):
        shutil.copytree(filename,new_name)
    else:
        shutil.copyfile(filename,new_name)
    result_files.append(new_name)        
    #write the report.html
    from genomicode import parselib
    from genomicode import htmllib
    def highlight(s):
        return htmllib.SPAN(s, style="background-color:yellow")
    def smaller(s):
        return htmllib.FONT(s, size=-1)
    try:
        lines = []
        w = lines.append
        w("<HTML>")
        title = "Heatmap Results"
        x = parselib.remove_all_tags(title)
        w(htmllib.HEAD(htmllib.TITLE(x)))
        w("<BODY>")
        w(htmllib.CENTER(htmllib.H1(title)))
        w(htmllib.P())
        w(htmllib.A("Methods",href="#methods_clustering"))
        w(htmllib.P())
        w(htmllib.A(htmllib.IMG(height=500,
            src=result_files[0]), href=result_files[0]))
        w(htmllib.P())
        name = 'Figure 1: In this heatmap, each row contains a signature and each column \
        contains a sample from your data set.'
        w(htmllib.B(name))
        
        w(htmllib.HR())
        w(htmllib.A("<methods_clustering>",name="methods_clustering"))
        w(htmllib.CENTER(htmllib.H2("Methods")))
        w(htmllib.H3("1.Heatmap File"))
        w('To generate this file, I ran the following analysis:')
        bie.plot_network_gv("network.png", network)
        w(htmllib.A(htmllib.IMG(height=500,
            src="network.png"), href="network.png"))
        w(htmllib.P())
        w('I used the following parameters:')
        w(htmllib.H3("1. Heatmap File"))
        rows = []
        x = htmllib.TR(
            htmllib.TH("Parameter", align="LEFT") +
            htmllib.TH("Value", align="LEFT") 
            )
        rows.append(x)
        
        for key in data_node.attributes.keys():
            x = htmllib.TR(
            htmllib.TD(key, align="LEFT") +
            htmllib.TD(data_node.attributes[key], align="LEFT") 
            )
            rows.append(x)
        w(htmllib.TABLE("\n".join(rows), border=1, cellpadding=3, cellspacing=0))
        w(htmllib.P())
        # Write out the footer.
        time_str = parselib.pretty_date(time.time())
        #hostname = pretty_hostname()
        w(htmllib.P())
        w(htmllib.HR())
        #w(htmllib.EM(
        #    "This analysis was run on %s on %s. \n" %
        #    (time_str, hostname)))
        w("</BODY>")
        w("</HTML>")
        x = "\n".join(lines) + "\n"
        open(outfile, 'w').write(x)
    except:
        raise
    
    new_parameters = parameters.copy()
    new_parameters['filename'] = os.path.split(outfile)[-1]
    out_node = bie.Data(rulebase.ReportFile,**new_parameters)
    return out_node

def name_outfile(data_node):
    filename = 'report.html' 
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(data_node,pipeline,parameters):
    identifier = data_node.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node = module_utils.get_identifier(network, module_id, data_nodes)
    return data_node
