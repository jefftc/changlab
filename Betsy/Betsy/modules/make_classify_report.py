#make_classify_report.py
from Betsy import config
from Betsy import module_utils
from Betsy import hash_method
import os
import shutil
from Betsy import protocol_utils
import imghdr
import time
import arrayio
import math
from genomicode import parselib
from genomicode import htmllib


def run(outfiles, parameters, pipelines):
    OUTPUTPATH = config.OUTPUTPATH
    inputid = module_utils.get_inputid(outfiles[0])
    folder_string = hash_method.hash_parameters(
            inputid, pipelines[0], **parameters[0])
    folder_name = 'classification_report_BETSYHASH1_' + folder_string
    result_folder = os.path.join(OUTPUTPATH, folder_name)
    if not os.path.exists(result_folder):
        os.mkdir(result_folder)
    result_files = []
    new_outfiles = []
    for j in range(len(outfiles)):
        filename = os.path.split(outfiles[j])[-1]
        folder = os.path.split(os.path.split(outfiles[j])[0])[-1]
        final_output = folder + '_' + filename
        result_files.append(final_output)
        result_file = os.path.join(result_folder, final_output)
        if not os.path.exists(result_file):
            if os.path.isdir(outfiles[j]):
                shutil.copytree(outfiles[j], result_file)
            else:
                shutil.copyfile(outfiles[j], result_file)

    #write the report.html
    
    def highlight(s):
        return htmllib.SPAN(s, style="background-color:yellow")
    def smaller(s):
        return htmllib.FONT(s, size=-1)
    
    cwd = os.getcwd()
    os.chdir(result_folder)
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
        M = arrayio.read(result_files[0])
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
        nfeature = int(parameters[6]['num_features'])
        if nfeature == 0:
            nfeature = 10
        M = arrayio.read(result_files[0])
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
        w(htmllib.CENTER(htmllib.H2("SVM Methods")))
        w(htmllib.H3("Prediction Result"))
        w('To generate this file, I ran the following analysis:')
        w(htmllib.P())
        for i in range(len(pipelines[1])):
            w('&nbsp&nbsp &nbsp&nbsp &nbsp&nbsp &nbsp&nbsp' +str(i+1)+'. '+pipelines[1][i])
            w(htmllib.P())
        w('I used the following parameters:')
        rows = []
        x = htmllib.TR(
            htmllib.TH("Parameter", align="LEFT") +
            htmllib.TH("Value", align="LEFT") 
            )
        rows.append(x)
        for key in parameters[1].keys():
            x = htmllib.TR(
            htmllib.TD(key, align="LEFT") +
            htmllib.TD(parameters[1][key], align="LEFT") 
            )
            rows.append(x)
        w(htmllib.TABLE("\n".join(rows), border=1, cellpadding=3, cellspacing=0))
        w(htmllib.P())
        w(htmllib.A("<methods_wv>",name="methods_wv"))
        w(htmllib.CENTER(htmllib.H2("Weighted Voting Methods")))
        w(htmllib.H3("Prediction Result"))
        w('To generate this file, I ran the following analysis:')
        w(htmllib.P())
        for i in range(len(pipelines[6])):
            w('&nbsp&nbsp &nbsp&nbsp &nbsp&nbsp &nbsp&nbsp' +str(i+1)+'. '+pipelines[6][i])
            w(htmllib.P())
        w('I used the following parameters:')
        rows = []
        x = htmllib.TR(
            htmllib.TH("Parameter", align="LEFT") +
            htmllib.TH("Value", align="LEFT") 
            )
        rows.append(x)
        for key in parameters[6].keys():
            x = htmllib.TR(
            htmllib.TD(key, align="LEFT") +
            htmllib.TD(parameters[6][key], align="LEFT") 
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
    finally:
        os.chdir(cwd)
    
    print 'Report:'+ os.path.join(result_folder,'report.html')

