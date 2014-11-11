#make_classify_report.py
import os
import shutil
import imghdr
import time
import arrayio
import math
from genomicode import parselib
from genomicode import htmllib
from Betsy import bie3
from Betsy import rulebase
from Betsy import module_utils
from Betsy import hash_method

def run(in_nodes, parameters, user_input, network,num_cores):
    outfile_folder=name_outfile(in_nodes,user_input)
    outfile = os.path.join(outfile_folder,'report.html')
    if not os.path.exists(outfile_folder):
        os.mkdir(outfile_folder)
    result_files = []
    for data_node in in_nodes:
        filename = data_node.identifier
        new_name = os.path.join(outfile_folder,os.path.split(filename)[-1])
        if os.path.isdir(filename):
                shutil.copytree(filename,new_name)
        else:
                shutil.copyfile(filename,new_name)
        result_files.append(os.path.split(new_name)[-1])   
    (data_node1, data_node2, data_node3, data_node4,data_node5, data_node6,
     data_node7,data_node8, data_node9, data_node10, data_node11) = in_nodes

    #write the report.html
    
    def highlight(s):
        return htmllib.SPAN(s, style="background-color:yellow")
    def smaller(s):
        return htmllib.FONT(s, size=-1)
    try:
        lines = []
        w = lines.append
        w("<HTML>")
        title = "Classification Results"
        x = parselib.remove_all_tags(title)
        w(htmllib.HEAD(htmllib.TITLE(x)))
        w("<BODY>")
        w(htmllib.CENTER(htmllib.H1(title)))
        #------------------------------------    
        w(htmllib.H3("SVM"))
        w(htmllib.P())
        w(htmllib.A("Methods",href="#methods_svm"))
        w(htmllib.P())
        #------------------------------------
        whole_row = []
        name = 'Table 1: Table of genes used in classification'
        w(htmllib.B(name))
        w(htmllib.P())
        M = arrayio.read(os.path.join(outfile_folder,result_files[0]))
        ids = M._row_order
        genes = M.row_names(ids[0])
        ncolumn = 3
        nrow = 8
        rows = []
        for i in range(min(nrow, len(genes)/ncolumn)):
            a= [] 
            for j in range(0, ncolumn):
                a.append('<td>' + genes[ncolumn * i + j] + '</td>')
            x = htmllib.TR("\n".join(a))
            rows.append(x)
        more_genes = 0
        if len(genes)>ncolumn*nrow:
            more_genes = len(genes)-ncolumn*nrow
        y = htmllib.TR(htmllib.TD(htmllib.TABLE("\n".join(rows), border=1, cellpadding=3, cellspacing=0),align='CENTER')+
                htmllib.TD(htmllib.A(htmllib.IMG(height=400,
                src=result_files[5]), href=result_files[5]),align='CENTER'))      
         #---------------------------------
        whole_row.append(y)
        y = htmllib.TR(htmllib.TD(htmllib.A(str(more_genes)+' more genes', result_files[0]),align='LEFT')+
                      htmllib.TD(htmllib.B('Figure 1: This figure shows the PCA plot of samples colored by prediction'),
                                 align='CENTER'))
        whole_row.append(y)
        x = htmllib.TR(
            htmllib.TD(htmllib.A(htmllib.IMG(height=400,
            src=result_files[4]), href=result_files[4]), align="CENTER") +
            htmllib.TD(htmllib.A(htmllib.IMG(height=400,
            src=result_files[2]), href=result_files[2]), align="CENTER") 
            )
        whole_row.append(x)
        x = htmllib.TR(
        htmllib.TH(htmllib.A("Figure 2. Loocv result on training data",result_files[3]), align="CENTER") +
        htmllib.TH(htmllib.A("Figure 3. Prediction result on test data",result_files[1]), align="CENTER") 
        )
        whole_row.append(x)
        w(htmllib.TABLE("\n".join(whole_row), border=None, cellpadding=3, cellspacing=0))
        w(htmllib.P())

        #------------------------------------    
        w(htmllib.H3("Weighted Voting"))
        w(htmllib.P())
        w(htmllib.A("Methods",href="#methods_wv"))
        w(htmllib.P())
        #------------------------------------
        whole_row = []
        name = 'Table 1: Table of genes used in classification'
        w(htmllib.B(name))
        w(htmllib.P())
        nfeature = 10
        if 'num_features_value' in user_input:
            nfeature = user_input['num_features_value']
        
        M = arrayio.read(os.path.join(outfile_folder,result_files[0]))
        ids = M._row_order
        genes = M.row_names(ids[0])[0:nfeature]
        nrow = min(8,int(math.ceil(float(len(genes))/ncolumn)))
        ncolumn = 3
        if len(genes)<nrow*ncolumn:
            genes.extend(['']*(nrow*ncolumn-len(genes)))
        rows = []
        for i in range(nrow):
            a= []
            for j in range(ncolumn):
                a.append('<td>'+genes[ncolumn*i+j]+'</td>')
            x = htmllib.TR("\n".join(a))
            rows.append(x)
        more_genes = 0
        if len(genes)>ncolumn*nrow:
            more_genes = len(genes)-ncolumn*nrow
            
        y = htmllib.TR(htmllib.TD(htmllib.TABLE("\n".join(rows), border=1, cellpadding=3, cellspacing=0),align='CENTER')+
                htmllib.TD(htmllib.A(htmllib.IMG(height=400,
                src=result_files[10]), href=result_files[10]),align='CENTER'))      
         #---------------------------------
        whole_row.append(y)
        y = htmllib.TR(htmllib.TD(htmllib.A(str(more_genes)+' more genes', result_files[0]),align='LEFT')+
                      htmllib.TD(htmllib.B('Figure 4: This figure shows the PCA plot of samples colored by prediction'),
                                 align='CENTER'))
        whole_row.append(y)
        x = htmllib.TR(
            htmllib.TD(htmllib.A(htmllib.IMG(height=400,
            src=result_files[9]), href=result_files[9]), align="CENTER") +
            htmllib.TD(htmllib.A(htmllib.IMG(height=400,
            src=result_files[7]), href=result_files[7]), align="CENTER") 
            )
        whole_row.append(x)
        x = htmllib.TR(
        htmllib.TH(htmllib.A("Figure 2. Loocv result on training data",result_files[8]), align="CENTER") +
        htmllib.TH(htmllib.A("Figure 3. Prediction result on test data",result_files[6]), align="CENTER") 
        )
        whole_row.append(x)
        w(htmllib.TABLE("\n".join(whole_row), border=None, cellpadding=3, cellspacing=0))
        w(htmllib.P())

       
        #--------------------------------
        
        w(htmllib.HR())
        w(htmllib.A("<methods_svm>",name="methods_svm"))
        w('To generate these files, I ran the following analysis:')
        bie3.plot_network_gv(os.path.join(outfile_folder,"network.png"), network)
        w(htmllib.P())
        w(htmllib.A(htmllib.IMG(height=500,
            src="network.png"), href="network.png"))
        w(htmllib.CENTER(htmllib.H2("SVM Methods")))
        w(htmllib.H3("Prediction Result"))
        
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
        w(htmllib.A("<methods_wv>",name="methods_wv"))
        w(htmllib.CENTER(htmllib.H2("Weighted Voting Methods")))
        w(htmllib.H3("Prediction Result"))
        w('I used the following parameters:')
        rows = []
        x = htmllib.TR(
            htmllib.TH("Parameter", align="LEFT") +
            htmllib.TH("Value", align="LEFT") 
            )
        rows.append(x)
        for key in data_node7.data.attributes.keys():
            x = htmllib.TR(
            htmllib.TD(key, align="LEFT") +
            htmllib.TD(data_node7.data.attributes[key], align="LEFT") 
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
    out_node = bie3.Data(rulebase.ClassifyReportFile,**parameters)
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
     data_node7,data_node8, data_node9, data_node10, data_node11) = in_nodes
    identifier = data_node1.identifier
    return module_utils.make_unique_hash(identifier,pipeline,parameters,user_input)

def find_antecedents(network, module_id,data_nodes,parameters,user_attributes):
    data_node1 = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='SignalFile')
    data_node2 = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='ClassifyFile',
                                             optional_key='classify_alg',
                                             optional_value='svm',second_key='loocv',
                                             second_value='no')
    data_node3 = module_utils.get_identifier(network, module_id, data_nodes,user_attributes,
                                            datatype='PredictionPlot',
                                             optional_key='classify_alg',
                                             optional_value='svm',second_key='loocv',
                                             second_value='no')
    data_node4 = module_utils.get_identifier(network, module_id, data_nodes,user_attributes,
                                            datatype='ClassifyFile',
                                             optional_key='classify_alg',
                                             optional_value='svm',second_key='loocv',
                                             second_value='yes')
    data_node5 = module_utils.get_identifier(network, module_id, data_nodes,user_attributes,
                                            datatype='PredictionPlot',
                                             optional_key='classify_alg',
                                             optional_value='svm',second_key='loocv',
                                             second_value='yes')
    data_node6 = module_utils.get_identifier(network, module_id, data_nodes,user_attributes,
                                           datatype='PredictionPCAPlot',
                                             optional_key='classify_alg',
                                             optional_value='svm')
    data_node7 = module_utils.get_identifier(network, module_id,
                                            data_nodes,user_attributes,datatype='ClassifyFile',
                                             optional_key='classify_alg',
                                             optional_value='weighted_voting',second_key='loocv',
                                             second_value='no')
    data_node8 = module_utils.get_identifier(network, module_id, data_nodes,user_attributes,
                                            datatype='PredictionPlot',
                                             optional_key='classify_alg',
                                             optional_value='weighted_voting',second_key='loocv',
                                             second_value='no')
    data_node9 = module_utils.get_identifier(network, module_id, data_nodes,user_attributes,
                                            datatype='ClassifyFile',
                                             optional_key='classify_alg',
                                             optional_value='weighted_voting',second_key='loocv',
                                             second_value='yes')
    data_node10 = module_utils.get_identifier(network, module_id, data_nodes,user_attributes,
                                            datatype='PredictionPlot',
                                             optional_key='classify_alg',
                                             optional_value='weighted_voting',second_key='loocv',
                                             second_value='yes')
    data_node11 = module_utils.get_identifier(network, module_id, data_nodes,user_attributes,
                                           datatype='PredictionPCAPlot',
                                             optional_key='classify_alg',
                                             optional_value='weighted_voting')
    
    return (data_node1, data_node2, data_node3, data_node4,data_node5, data_node6,
            data_node7,data_node8, data_node9,data_node10, data_node11)

