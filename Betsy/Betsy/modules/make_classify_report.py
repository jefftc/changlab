#make_classify_report.py

import config
import module_utils
import hash_method
import os
import shutil
import protocol_utils
import imghdr
import time
import arrayio
from genomicode import parselib
from genomicode import htmllib
def run(outfiles,parameters,pipelines):
    OUTPUTPATH = config.OUTPUTPATH
    inputid = module_utils.get_inputid(outfiles[0])
    folder_string = hash_method.hash_parameters(
            inputid,pipelines[0],**parameters[0])
    folder_name = 'classification_report_BETSYHASH1_'+folder_string
    result_folder = os.path.join(OUTPUTPATH,folder_name)
    if not os.path.exists(result_folder):
        os.mkdir(result_folder)
    result_files = []
    new_outfiles = []
    for j in range(len(outfiles)):
        filename = os.path.split(outfiles[j])[-1]
        folder = os.path.split(os.path.split(outfiles[j])[0])[-1]
        final_output = folder+'_'+ filename
        result_files.append(final_output)
        result_file = os.path.join(result_folder,final_output)
        if not os.path.exists(result_file):
            if os.path.isdir(outfiles[j]):
                shutil.copytree(outfiles[j],result_file)
            else:
                shutil.copyfile(outfiles[j],result_file)

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
        #w('I generated a file that contains the normalized gene expression values')
        #w(htmllib.P())
        #w(htmllib.A(result_files[0], result_files[0]))
        name = 'Table 1: Table of genes used in classification'
        w(htmllib.B(name))
        w(htmllib.P())
        M = arrayio.read(result_files[0])
        ids = M._row_order
        genes = M.row_names(ids[0])
        rows = []
        a = htmllib.TH('Genes',align='LEFT')
        x = htmllib.TR(a)
        rows.append(x)
        for i in range(min(50,len(genes))):
            a=htmllib.TD(genes[i],align="LEFT")
            x = htmllib.TR(a)
            rows.append(x)
        #rows = write_table(header,[genes],50)
        w(htmllib.TABLE("\n".join(rows), border=1, cellpadding=3, cellspacing=0))
        w(htmllib.P())
        if len(genes)>51:
            more_genes = len(genes)-51
            w(htmllib.A(str(more_genes)+' more genes', result_files[0]))
            w(htmllib.P())
            
        w(htmllib.H3("SVM"))
        w(htmllib.P())
        w(htmllib.A("Methods",href="#methods_svm"))
        w(htmllib.P())
         #---------------------------------
        name = 'Table 2: Table of prediction result of LOOCV on training set'
        w(htmllib.B(name))
        f = file(result_files[1],'rU')
        text = f.readlines()
        f.close()
        header = text[0].split('\t')
        data=[i.split('\t') for i in text[1:]]
        rows = write_table(header,data,50)
        w(htmllib.TABLE("\n".join(rows), border=1, cellpadding=3, cellspacing=0))
        w(htmllib.P())
        if len(text)>51:
            more_genes = len(text)-51
            w(htmllib.A(str(more_genes)+' more genes', result_files[1]))
            w(htmllib.P())
        
        #---------------------------------
        name = 'Table 3: Table of prediction result of on test set'
        w(htmllib.B(name))
        f = file(result_files[2],'rU')
        text = f.readlines()
        f.close()
        header = text[0].split('\t')
        data=[i.split('\t') for i in text[1:]]
        rows = write_table(header,data,50)
        w(htmllib.TABLE("\n".join(rows), border=1, cellpadding=3, cellspacing=0))
        w(htmllib.P())
        if len(text)>51:
            more_genes = len(text)-51
            w(htmllib.A(str(more_genes)+' more genes', result_files[2]))
            w(htmllib.P())
        
        #---------------------------------
        name = 'Figure 1: This figure shows the PCA plot of samples colored by prediction'
        w(htmllib.B(name))
        w(htmllib.P())
        w(htmllib.A(htmllib.IMG(height=500,
                src=result_files[3]), href=result_files[3]))
        w(htmllib.P())
        #---------------------------------
        w(htmllib.H3("Weighted Voting"))
        w(htmllib.A("Methods",href="#methods_wv"))
        w(htmllib.P())
        name = 'Table 4: Table of prediction result of LOOCV on training set'
        w(htmllib.B(name))
        f = file(result_files[4],'rU')
        text = f.readlines()
        f.close()
        header = text[0].split('\t')
        data=[i.split('\t') for i in text[1:]]
        rows = write_table(header,data,50)
        w(htmllib.TABLE("\n".join(rows), border=1, cellpadding=3, cellspacing=0))
        w(htmllib.P())
        if len(text)>51:
            more_genes = len(text)-51
            w(htmllib.A(str(more_genes)+' more genes', result_files[4]))
            w(htmllib.P())
        
        #---------------------------------
        name = 'Table 3: Table of prediction result of on test set'
        w(htmllib.B(name))
        f = file(result_files[5],'rU')
        text = f.readlines()
        f.close()
        header = text[0].split('\t')
        data=[i.split('\t') for i in text[1:]]
        rows = write_table(header,data,50)
        w(htmllib.TABLE("\n".join(rows), border=1, cellpadding=3, cellspacing=0))
        w(htmllib.P())
        if len(text)>51:
            more_genes = len(text)-51
            w(htmllib.A(str(more_genes)+' more genes', result_files[5]))
            w(htmllib.P())
        
        #---------------------------------
        name = 'Figure 2: This figure shows the PCA plot of samples colored by prediction'
        w(htmllib.B(name))
        w(htmllib.P())
        w(htmllib.A(htmllib.IMG(height=500,
                src=result_files[6]), href=result_files[6]))
        w(htmllib.P())
        #--------------------------------
        
        w(htmllib.HR())
        w(htmllib.A("<methods_svm>",name="methods_svm"))
        w(htmllib.CENTER(htmllib.H2("SVM Methods")))
        w(htmllib.H3("Prediction Result"))
        w('To generate this file, I ran the following analysis:')
        w(htmllib.P())
        for i in range(len(pipelines[2])):
            w('&nbsp&nbsp &nbsp&nbsp &nbsp&nbsp &nbsp&nbsp' +str(i+1)+'. '+pipelines[2][i])
            w(htmllib.P())
        w('I used the following parameters:')
        rows = []
        x = htmllib.TR(
            htmllib.TH("Parameter", align="LEFT") +
            htmllib.TH("Value", align="LEFT") 
            )
        rows.append(x)
        for key in parameters[2].keys():
            x = htmllib.TR(
            htmllib.TD(key, align="LEFT") +
            htmllib.TD(parameters[2][key], align="LEFT") 
            )
            rows.append(x)
        w(htmllib.TABLE("\n".join(rows), border=1, cellpadding=3, cellspacing=0))
        w(htmllib.P())
        w(htmllib.A("<methods_wv>",name="methods_wv"))
        w(htmllib.CENTER(htmllib.H2("Weighted Voting Methods")))
        w(htmllib.H3("Prediction Result"))
        w('To generate this file, I ran the following analysis:')
        w(htmllib.P())
        for i in range(len(pipelines[5])):
            w('&nbsp&nbsp &nbsp&nbsp &nbsp&nbsp &nbsp&nbsp' +str(i+1)+'. '+pipelines[5][i])
            w(htmllib.P())
        w('I used the following parameters:')
        rows = []
        x = htmllib.TR(
            htmllib.TH("Parameter", align="LEFT") +
            htmllib.TH("Value", align="LEFT") 
            )
        rows.append(x)
        for key in parameters[5].keys():
            x = htmllib.TR(
            htmllib.TD(key, align="LEFT") +
            htmllib.TD(parameters[5][key], align="LEFT") 
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
