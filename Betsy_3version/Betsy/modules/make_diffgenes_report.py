#make_diffgenes_report.py
import os
import shutil
import imghdr
import time
from genomicode import parselib,htmllib
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils
from Betsy import hash_method

def run(in_nodes,parameters,user_input, network):
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
    data_node1, data_node2, data_node3, data_node4,data_node5 = in_nodes
    #write the report.html
    
    def highlight(s):
        return htmllib.SPAN(s, style="background-color:yellow")
    def smaller(s):
        return htmllib.FONT(s, size=-1)
    try:
        lines = []
        w = lines.append
        w("<HTML>")
        title = "Differential Expressed Gene Analysis Result"
        x = parselib.remove_all_tags(title)
        w(htmllib.HEAD(htmllib.TITLE(x)))
        w("<BODY>")
        w(htmllib.CENTER(htmllib.H1(title)))
        w('I generated a file that contains the t-test results')
        w(htmllib.P())
        w(htmllib.A(result_files[0], result_files[0]))
        w(htmllib.P())
        w(htmllib.A("Methods",href="#methods_diffgenes"))
        w(htmllib.P())
        w('I generated a file that contains the sam results')
        w(htmllib.P())
        sam_result = os.path.join(result_files[1],'sam_result.txt')
        sam_fig = os.path.join(result_files[1],'sam_plot.png')
        w(htmllib.A(sam_result, sam_result))
        w(htmllib.P())
        #---------------------------------
        name = 'Table 1: Table of significant genes p<0.05 sorted in order of significance'
        w(htmllib.B(name))
        f = file(result_files[0],'rU')
        text = f.readlines()
        f.close()
        header = text[0].split('\t')
        data=[i.split('\t') for i in text[1:]]
        rows = write_table(header,data,50)
        w(htmllib.TABLE("\n".join(rows), border=1, cellpadding=3, cellspacing=0))
        w(htmllib.P())
        if len(text)>51:
            more_genes = len(text)-51
            w(htmllib.A(str(more_genes)+' more genes', result_files[0]))
            w(htmllib.P())
        #-----------------------------------
        w(htmllib.A(htmllib.IMG(height=500,
                src=result_files[2]), href=result_files[2]))
        w(htmllib.P())
        name = 'Figure 1: Heatmap of significant genes'
        w(htmllib.B(name))
        w(htmllib.P())
        #-----------------------------------
        name = 'Table 2: Table of significant annotations'
        w(htmllib.B(name))
        w(htmllib.P())
        f = file(result_files[3],'rU')
        text = f.readlines()
        f.close()
        index = [0,1,2,3,4,5,6,7,9,10]
        header = text[0].split('\t')
        header = [header[i] for i in index]
        data=[i.split('\t') for i in text[1:]]
        data=[[item[i] for i in index] for item in data]
        import math
        p_list=[(data[i][8],i) for i in range(len(data)) if data[i][8]>-math.log(0.05)]
        p_list.sort(reverse=True)
        sort_index = [i[1] for i in p_list]
        data = [data[i] for i in sort_index]
        rows = write_table(header,data,10)
        w(htmllib.TABLE("\n".join(rows), border=1, cellpadding=3, cellspacing=0))
        w(htmllib.P())
        if len(data)>10:
            more_genes = len(data)-10
            w(htmllib.A(str(more_genes)+' more annotations', result_files[3]))
            w(htmllib.P())
        #-----------------------------------
        w(htmllib.A(htmllib.IMG(height=500,
                src=sam_fig), href=sam_fig))
        w(htmllib.P())
        name = 'Figure 2: SAM plot'
        w(htmllib.B(name))
        #-----------------------------------
        w(htmllib.P())
        w('The full result of Gather is in')
        w(htmllib.P())
        w(htmllib.A(result_files[3], result_files[3]))
        w(htmllib.P())
        #-----------------------------------
        w('The result of GSEA is in')
        w(htmllib.P())
        w(htmllib.A(result_files[4], result_files[4]))
        w(htmllib.P())
        w(htmllib.HR())
        #-----------------------------------
        w(htmllib.A("<methods_diffgenes>",name="methods_diffgenes"))
        w(htmllib.CENTER(htmllib.H2("Methods")))
        w(htmllib.H3("1.T-test"))
        w('To generate this file, I ran the following analysis:')
        bie.plot_network_gv("network.png", network)
        w(htmllib.P())
        w(htmllib.A(htmllib.IMG(height=500,
            src="network.png"), href="network.png"))
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
        w("</BODY>")
        w("</HTML>")
        x = "\n".join(lines) + "\n"
        open('report.html', 'w').write(x)
    except:
        raise 
    out_node = bie3.Data(rulebase.ReportFile,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object
    
    
def write_table(header,data,N):
    rows = []
    a = ''
    for i in header:
        a= a + htmllib.TH(i,align='LEFT')
    x = htmllib.TR(a)
    rows.append(x)
    for i in range(min(N,len(data))):
        a=''
        for j in range(len(data[0])):
             a=a+htmllib.TD(data[i][j],align="LEFT")
        x = htmllib.TR(a)
        rows.append(x)
    return rows

def name_outfile(in_nodes,user_input):
    filename = 'report.html' 
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters,user_input):
    data_node1, data_node2, data_node3, data_node4,data_node5 = in_nodes
    identifier = data_node1.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)

def find_antecedents(network, module_id,data_nodes,parameters):
    data_node1 = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='DiffExprFile',
                                             optional_key='diff_expr',
                                             optional_value='t_test')
    data_node2 = module_utils.get_identifier(network, module_id,
                                            data_nodes,datatype='DiffExprFile',
                                             optional_key='diff_expr',
                                             optional_value='sam')
    data_node3 = module_utils.get_identifier(network, module_id, data_nodes,
                                           datatype='Heatmap')
    data_node4 = module_utils.get_identifier(network, module_id, data_nodes,
                                           datatype='GatherFile')
    data_node5 = module_utils.get_identifier(network, module_id, data_nodes,
                                           datatype='GseaFile')
    return data_node1, data_node2, data_node3, data_node4,data_node5
