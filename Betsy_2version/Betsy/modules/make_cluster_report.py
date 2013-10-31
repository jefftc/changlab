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

def run(in_nodes,parameters, network):
    outfile=name_outfile(in_nodes)
    result_files = []
    for data_node in in_nodes:
        filename = data_node.attributes['filename']
        new_name = os.path.split(filename)[-1]
        if os.path.isdir(filename):
                shutil.copytree(filename,new_name)
        else:
                shutil.copyfile(filename,new_name)
        result_files.append(new_name)        
    data_node1,data_node2 = in_nodes
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
        title = "Clustering Results"
        x = parselib.remove_all_tags(title)
        w(htmllib.HEAD(htmllib.TITLE(x)))
        w("<BODY>")
        w(htmllib.CENTER(htmllib.H1(title)))
        w('I generated a file that contains the expression values of the data set after clustering')
        w(htmllib.P())
        w(htmllib.A(result_files[0], result_files[0]))
        w(htmllib.P())
        w(htmllib.A("Methods",href="#methods_clustering"))
        w(htmllib.P())
        w(htmllib.A(htmllib.IMG(height=500,
            src=result_files[1]), href=result_files[1]))
        w(htmllib.P())
        name = 'Figure 1: In this heatmap, each row contains a signature and each column \
        contains a sample from your data set.'
        w(htmllib.B(name))
        
        w(htmllib.HR())
        w(htmllib.A("<methods_clustering>",name="methods_clustering"))
        w(htmllib.CENTER(htmllib.H2("Methods")))
        w('To generate this file, I ran the following analysis:')
        bie._plot_network_gv("network.png", network)
##        w(htmllib.P())
##        for i in range(len(pipelines[0])):
##            w('&nbsp&nbsp &nbsp&nbsp &nbsp&nbsp &nbsp&nbsp' +str(i+1)+'. '+pipelines[0][i])
##            w(htmllib.P())
        w(htmllib.A(htmllib.IMG(height=500,
            src="network.png"), href="network.png"))
        w(htmllib.P())
        name = 'Figure 1: In this heatmap, each row contains a signature and each column \
        contains a sample from your data set.'
        w('I used the following parameters:')
        w(htmllib.H3("1.Cluster File"))
        rows = []
        x = htmllib.TR(
            htmllib.TH("Parameter", align="LEFT") +
            htmllib.TH("Value", align="LEFT") 
            )
        rows.append(x)
        data_node1.attributes
        for key in data_node1.attributes.keys():
            x = htmllib.TR(
            htmllib.TD(key, align="LEFT") +
            htmllib.TD(data_node1.attributes[key], align="LEFT") 
            )
            rows.append(x)
        w(htmllib.TABLE("\n".join(rows), border=1, cellpadding=3, cellspacing=0))
        w(htmllib.P())
        w(htmllib.H3("2.Cluster heatmap"))
        rows = []
        x = htmllib.TR(
            htmllib.TH("Parameter", align="LEFT") +
            htmllib.TH("Value", align="LEFT") 
            )
        rows.append(x)
        data_node1.attributes
        for key in data_node2.attributes.keys():
            x = htmllib.TR(
            htmllib.TD(key, align="LEFT") +
            htmllib.TD(data_node2.attributes[key], align="LEFT") 
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

def name_outfile(in_nodes):
    filename = 'report.html' 
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters):
    data_node1,data_node2 = in_nodes
    identifier = data_node1.attributes['filename']
    return module_utils.make_unique_hash(identifier,pipeline,parameters)

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node1 = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='ClusterFile')
    data_node2 = module_utils.get_identifier(network, module_id, data_nodes,
                                           datatype='Heatmap')
    return data_node1, data_node2