#make_geneset_report.py

from Betsy import module_utils
from Betsy import hash_method
import os
import shutil
import imghdr
import time
from genomicode import parselib
from genomicode import htmllib
from Betsy import bie3
from Betsy import rulebase

def run(in_nodes, parameters, user_input,network):
    outfile=name_outfile(in_nodes,user_input)
    result_files = []
    for data_node in in_nodes:
        filename = data_node.identifier
        new_name = os.path.split(filename)[-1]
        if os.path.isdir(filename):
                shutil.copytree(filename,new_name)
        else:
                shutil.copyfile(filename,new_name)
        result_files.append(new_name)        
    data_node1,data_node2 = in_nodes    #write the report.html
    
    def highlight(s):
        return htmllib.SPAN(s, style="background-color:yellow")
    def smaller(s):
        return htmllib.FONT(s, size=-1)
    
    
    try:
        lines = []
        w = lines.append
        w("<HTML>")
        title = "Geneset Analysis Results"
        x = parselib.remove_all_tags(title)
        w(htmllib.HEAD(htmllib.TITLE(x)))
        w("<BODY>")
        w(htmllib.CENTER(htmllib.H1(title)))
        w('I generated a file that contains the analysis result of the geneset')
        w(htmllib.P())
        w(htmllib.A(result_files[0], result_files[0]))
        w(htmllib.P())
        w(htmllib.A("Methods",href="#methods"))
        w(htmllib.P())
        filenames = os.listdir(result_files[1])
        c=0
        for filename in filenames:
            c=c+1
            w(htmllib.A(htmllib.IMG(height=500,
                src=os.path.join(result_files[1],filename)), href=os.path.join(result_files[1],filename)))
            w(htmllib.P())
            name = 'Figure '+ str(c) + ': Geneset Plot.'
            w(htmllib.B(name))
        
        w(htmllib.HR())
        w(htmllib.A("<methods>",name="methods"))
        w(htmllib.CENTER(htmllib.H2("Methods")))
        w(htmllib.H3("1.Result File"))
        w('To generate this file, I ran the following analysis:')
        bie.plot_network_gv("network.png", network)
        w(htmllib.A(htmllib.IMG(height=500,
            src="network.png"), href="network.png"))
        w(htmllib.P())
        
        w('I used the following parameters:')
        rows = []
        x = htmllib.TR(
            htmllib.TH("Parameter", align="LEFT") +
            htmllib.TH("Value", align="LEFT") 
            )
        rows.append(x)
        for key in data_node1.attributes.keys():
            x = htmllib.TR(
            htmllib.TD(key, align="LEFT") +
            htmllib.TD(data_node1.attributes[key], align="LEFT") 
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
        open('report.html', 'w').write(x)
    except:
        raise 
    out_node = bie3.Data(rulebase.ReportFile,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object
    
def name_outfile(in_nodes,user_input):
    filename = 'report.html' 
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters,user_input):
    data_node1,data_node2 = in_nodes
    identifier = data_node1.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)

def find_antecedents(network, module_id,data_nodes):
    data_node1 = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='GenesetAnalysis')
    data_node2 = module_utils.get_identifier(network, module_id, data_nodes,
                                           datatype='GenesetPlot')
    return data_node1, data_node2
