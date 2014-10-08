#make_batch_effect_report.py
import os
import shutil
import imghdr
import time
from genomicode import parselib,htmllib
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils
from Betsy import hash_method


def run(in_nodes, parameters, user_input, network):
    outfile_folder=name_outfile(in_nodes,user_input)
    outfile = os.path.join(outfile_folder,'report.html')
    if not os.path.exists(outfile_folder):
        os.mkdir(outfile_folder)
    result_files = []
    for data_node in in_nodes:
        filename = data_node.identifier
        folder = os.path.split(os.path.split(filename)[0])[-1]
        new_name = os.path.join(outfile_folder,folder+'_'+os.path.split(filename)[-1])
        if os.path.isdir(filename):
                shutil.copytree(filename,new_name)
        else:
                shutil.copyfile(filename,new_name)
        result_files.append(os.path.split(new_name)[-1])   
    (data_node1, data_node2, data_node3, data_node4,data_node5, data_node6,
     data_node7,data_node8, data_node9, data_node10, data_node11,date_node12) = in_nodes

    #write the report.html
    
    def highlight(s):
        return htmllib.SPAN(s, style="background-color:yellow")
    def smaller(s):
        return htmllib.FONT(s, size=-1)
    try:
        lines = []
        w = lines.append
        w("<HTML>")
        title = "Batch Effect Remove Results"
        x = parselib.remove_all_tags(title)
        w(htmllib.HEAD(htmllib.TITLE(x)))
        w("<BODY>")
        w(htmllib.CENTER(htmllib.H1(title)))
        w(htmllib.P())
        w(htmllib.A("Methods",href="#methods_batch_effect"))
        w(htmllib.P())
        w('I generated a file that contains the gene expression values of your dataset before batch effect remove')
        w(htmllib.P())
        w(htmllib.A(result_files[0], result_files[0]))
        w(htmllib.P())
        rows=[]
        x = htmllib.A(htmllib.IMG(height=500,
            src=result_files[1]), href=result_files[1])
        rows.append(x)
        
        w(htmllib.TABLE("\n".join(rows), border=None, cellpadding=3, cellspacing=0))
        w(htmllib.P())

        w(htmllib.P())
        name = 'Figure 1: This pca plot shows the similarities among your samples before batch effect remove'
        w(htmllib.B(name))
        w(htmllib.P())
    #----------------------------------------------------------------------------
        w(htmllib.P())
        w('I generated a file that contains the gene expression values of your dataset after quantile')
        w(htmllib.P())
        w(htmllib.A(result_files[2], result_files[2]))
        w(htmllib.P())
        rows=[]
        x = htmllib.A(htmllib.IMG(height=500,
            src=result_files[3]), href=result_files[3])
        rows.append(x)
        
        w(htmllib.TABLE("\n".join(rows), border=None, cellpadding=3, cellspacing=0))
        w(htmllib.P())

        w(htmllib.P())
        name = 'Figure 2: This pca plot shows the similarities among your samples after quantile'
        w(htmllib.B(name))
        w(htmllib.P())
        #----------------------------------------------------------------------------
        w(htmllib.P())
        w('I generated a file that contains the gene expression values of your dataset after quantile and bfrm')
        w(htmllib.P())
        w(htmllib.A(result_files[4], result_files[4]))
        w(htmllib.P())
        rows=[]
        x = htmllib.A(htmllib.IMG(height=500,
            src=result_files[5]), href=result_files[5])
        rows.append(x)
        
        w(htmllib.TABLE("\n".join(rows), border=None, cellpadding=3, cellspacing=0))
        w(htmllib.P())

        w(htmllib.P())
        name = 'Figure 3: This pca plot shows the similarities among your samples after quantile and bfrm'
        w(htmllib.B(name))
        w(htmllib.P())
        #----------------------------------------------------------------------------
        
        w(htmllib.P())
        w('I generated a file that contains the gene expression values of your dataset after quantile and combat')
        w(htmllib.P())
        w(htmllib.A(result_files[6], result_files[6]))
        w(htmllib.P())
        rows=[]
        x = htmllib.A(htmllib.IMG(height=500,
            src=result_files[7]), href=result_files[7])
        rows.append(x)
        
        w(htmllib.TABLE("\n".join(rows), border=None, cellpadding=3, cellspacing=0))
        w(htmllib.P())

        w(htmllib.P())
        name = 'Figure 4: This pca plot shows the similarities among your samples after quantile and combat'
        w(htmllib.B(name))
        w(htmllib.P())
        #----------------------------------------------------------------------------
        w(htmllib.P())
        w('I generated a file that contains the gene expression values of your dataset after quantile and dwd')
        w(htmllib.P())
        w(htmllib.A(result_files[8], result_files[8]))
        w(htmllib.P())
        rows=[]
        x = htmllib.A(htmllib.IMG(height=500,
            src=result_files[9]), href=result_files[9])
        rows.append(x)
        
        w(htmllib.TABLE("\n".join(rows), border=None, cellpadding=3, cellspacing=0))
        w(htmllib.P())

        w(htmllib.P())
        name = 'Figure 5: This pca plot shows the similarities among your samples after quantile and dwd'
        w(htmllib.B(name))
        w(htmllib.P())
        #----------------------------------------------------------------------------
        w(htmllib.P())
        w('I generated a file that contains the gene expression values of your dataset after quantile and shiftscale')
        w(htmllib.P())
        w(htmllib.A(result_files[10], result_files[10]))
        w(htmllib.P())
        rows=[]
        x = htmllib.A(htmllib.IMG(height=500,
            src=result_files[11]), href=result_files[11])
        rows.append(x)
        
        w(htmllib.TABLE("\n".join(rows), border=None, cellpadding=3, cellspacing=0))
        w(htmllib.P())

        w(htmllib.P())
        name = 'Figure 5: This pca plot shows the similarities among your samples after quantile and shiftscale'
        w(htmllib.B(name))
        w(htmllib.P())

    
        #--------------------------------
        
        w(htmllib.HR())
        w(htmllib.A("<methods_batch_effect>",name="methods"))
        w('To generate these files, I ran the following analysis:')
        bie3.plot_network_gv(os.path.join(outfile_folder,"network.png"), network)
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
        for key in data_node2.data.attributes.keys():
            x = htmllib.TR(
            htmllib.TD(key, align="LEFT") +
            htmllib.TD(data_node2.data.attributes[key], align="LEFT") 
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
    out_node = bie3.Data(rulebase.BatchEffectReportFile,**parameters)
    out_object = module_utils.DataObject(out_node,outfile)
    return out_object
    


def name_outfile(in_nodes,user_input):
    filename = 'report' 
    outfile = os.path.join(os.getcwd(), filename)
    return outfile

    
def get_out_attributes(parameters,in_nodes):
    return parameters

def make_unique_hash(in_nodes,pipeline,parameters,user_input):
    (data_node1, data_node2, data_node3, data_node4,data_node5, data_node6,
     data_node7,data_node8, data_node9, data_node10, data_node11,data_node12) = in_nodes
    identifier = data_node1.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node1 = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,
                                             datatype='SignalFile',
                                             optional_key='quantile_norm',
                                             optional_value='no')
    data_node2 = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,
                                             datatype='PcaPlot',
                                             optional_key='quantile_norm',
                                             optional_value='no')
    param = {'shiftscale_norm':'no',
             'dwd_norm':'no',
             'combat_norm':'no'}
    data_node3 = module_utils.get_identifier(network, module_id, data_nodes,user_attributes,
                                            datatype='SignalFile',
                                             optional_key='quantile_norm',
                                             optional_value='yes',second_key='bfrm_norm',
                                             second_value='no',**param) 
    data_node4 = module_utils.get_identifier(network, module_id, data_nodes,user_attributes,
                                            datatype='PcaPlot',
                                             optional_key='quantile_norm',
                                             optional_value='yes',second_key='bfrm_norm',
                                             second_value='no',**param) 
    data_node5 = module_utils.get_identifier(network, module_id, data_nodes,user_attributes,
                                            datatype='SignalFile',
                                             optional_key='quantile_norm',
                                             optional_value='yes',second_key='bfrm_norm',
                                             second_value='yes')
    data_node6 = module_utils.get_identifier(network, module_id, data_nodes,user_attributes,
                                             datatype='PcaPlot',
                                             optional_key='quantile_norm',
                                             optional_value='yes',second_key='bfrm_norm',
                                             second_value='yes')
    data_node7 = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='SignalFile',
                                             optional_key='quantile_norm',
                                             optional_value='yes',second_key='combat_norm',
                                             second_value='yes')
    data_node8 = module_utils.get_identifier(network, module_id, data_nodes,user_attributes,
                                            datatype='PcaPlot',
                                             optional_key='quantile_norm',
                                             optional_value='yes',second_key='combat_norm',
                                             second_value='yes')
    data_node9 = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='SignalFile',
                                             optional_key='quantile_norm',
                                             optional_value='yes',second_key='dwd_norm',
                                             second_value='yes')
    data_node10 = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='PcaPlot',
                                             optional_key='quantile_norm',
                                             optional_value='yes',second_key='dwd_norm',
                                             second_value='yes')
    data_node11 = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='SignalFile',
                                             optional_key='quantile_norm',
                                             optional_value='yes',second_key='shiftscale_norm',
                                             second_value='yes')
    data_node12 = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='PcaPlot',
                                             optional_key='quantile_norm',
                                             optional_value='yes',second_key='shiftscale_norm',
                                             second_value='yes')
    return (data_node1, data_node2, data_node3, data_node4,data_node5, data_node6,
            data_node7,data_node8, data_node9,data_node10, data_node11,data_node12)
