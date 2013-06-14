#make_batchnorm_report.py
from Betsy import config
from Betsy import module_utils
from Betsy import hash_method
import os
import shutil
from Betsy import protocol_utils
import imghdr
import time
from time import strftime,localtime

def run(outfiles,parameters,pipelines,user,jobname):
    starttime = strftime(module_utils.FMT, localtime())
    OUTPUTPATH = config.OUTPUTPATH
    inputid = module_utils.get_inputid(outfiles[0])
    folder_string = hash_method.hash_parameters(
            inputid,pipelines[0],**parameters[0])
    folder_name = 'batch_effect_remove_report_BETSYHASH1_'+folder_string
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
    from genomicode import parselib
    from genomicode import htmllib
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
        title = "Batch Effect Remove Results"
        x = parselib.remove_all_tags(title)
        w(htmllib.HEAD(htmllib.TITLE(x)))
        w("<BODY>")
        w(htmllib.CENTER(htmllib.H1(title)))
        w(htmllib.P())
        w(htmllib.A("Methods",href="#methods_normalization"))
        w(htmllib.P())
        w('I generated a file that contains the gene expression values of your dataset after quantile')
        w(htmllib.P())
        w(htmllib.A(result_files[0], result_files[0]))
        w(htmllib.P())
        rows=[]
        x = htmllib.TR(
            htmllib.TD(htmllib.A(htmllib.IMG(height=500,
            src=result_files[2]), href=result_files[2]), align="CENTER") +
            htmllib.TD(htmllib.A(htmllib.IMG(height=500,
            src=result_files[1]), href=result_files[1]), align="CENTER") 
            )
        rows.append(x)
        x = htmllib.TR(
        htmllib.TH("Before", align="CENTER") +
        htmllib.TH("After", align="CENTER") 
        )
        rows.append(x)
        w(htmllib.TABLE("\n".join(rows), border=None, cellpadding=3, cellspacing=0))
        w(htmllib.P())

        w(htmllib.P())
        name = 'Figure 1: This pca plot shows the similarities among your samples after quantile'
        w(htmllib.B(name))
        w(htmllib.P())
    #----------------------------------------------------------------------------
        w(htmllib.P())
        w('I generated a file that contains the gene expression values of your dataset after bfrm')
        w(htmllib.P())
        w(htmllib.A(result_files[3], result_files[3]))
        w(htmllib.P())
        rows=[]
        x = htmllib.TR(
            htmllib.TD(htmllib.A(htmllib.IMG(height=500,
            src=result_files[5]), href=result_files[5]), align="CENTER") +
            htmllib.TD(htmllib.A(htmllib.IMG(height=500,
            src=result_files[4]), href=result_files[4]), align="CENTER") 
            )
        rows.append(x)
        x = htmllib.TR(
        htmllib.TH("Before", align="CENTER") +
        htmllib.TH("After", align="CENTER") 
        )
        rows.append(x)
        w(htmllib.TABLE("\n".join(rows), border=None, cellpadding=3, cellspacing=0))
        w(htmllib.P())

        w(htmllib.P())
        name = 'Figure 2: This pca plot shows the similarities among your samples after bfrm'
        w(htmllib.B(name))
        w(htmllib.P())
        #----------------------------------------------------------------------------
        w(htmllib.P())
        w('I generated a file that contains the gene expression values of your dataset after quantile and combat')
        w(htmllib.P())
        w(htmllib.A(result_files[6], result_files[6]))
        w(htmllib.P())
        rows=[]
        x = htmllib.TR(
            htmllib.TD(htmllib.A(htmllib.IMG(height=500,
            src=result_files[8]), href=result_files[8]), align="CENTER") +
            htmllib.TD(htmllib.A(htmllib.IMG(height=500,
            src=result_files[7]), href=result_files[7]), align="CENTER") 
            )
        rows.append(x)
        x = htmllib.TR(
        htmllib.TH("Before", align="CENTER") +
        htmllib.TH("After", align="CENTER") 
        )
        rows.append(x)
        w(htmllib.TABLE("\n".join(rows), border=None, cellpadding=3, cellspacing=0))
        w(htmllib.P())

        w(htmllib.P())
        name = 'Figure 3: This pca plot shows the similarities among your samples after quantile and combat'
        w(htmllib.B(name))
        w(htmllib.P())
        #----------------------------------------------------------------------------
        
        w(htmllib.P())
        w('I generated a file that contains the gene expression values of your dataset after quantile and dwd')
        w(htmllib.P())
        w(htmllib.A(result_files[9], result_files[9]))
        w(htmllib.P())
        rows=[]
        x = htmllib.TR(
            htmllib.TD(htmllib.A(htmllib.IMG(height=500,
            src=result_files[11]), href=result_files[11]), align="CENTER") +
            htmllib.TD(htmllib.A(htmllib.IMG(height=500,
            src=result_files[10]), href=result_files[10]), align="CENTER") 
            )
        rows.append(x)
        x = htmllib.TR(
        htmllib.TH("Before", align="CENTER") +
        htmllib.TH("After", align="CENTER") 
        )
        rows.append(x)
        w(htmllib.TABLE("\n".join(rows), border=None, cellpadding=3, cellspacing=0))
        w(htmllib.P())

        w(htmllib.P())
        name = 'Figure 4: This pca plot shows the similarities among your samples after quantile and dwd'
        w(htmllib.B(name))
        w(htmllib.P())
        #----------------------------------------------------------------------------
        w(htmllib.P())
        w('I generated a file that contains the gene expression values of your dataset after quantile and shiftscale')
        w(htmllib.P())
        w(htmllib.A(result_files[12], result_files[12]))
        w(htmllib.P())
        rows=[]
        x = htmllib.TR(
            htmllib.TD(htmllib.A(htmllib.IMG(height=500,
            src=result_files[14]), href=result_files[14]), align="CENTER") +
            htmllib.TD(htmllib.A(htmllib.IMG(height=500,
            src=result_files[13]), href=result_files[13]), align="CENTER") 
            )
        rows.append(x)
        x = htmllib.TR(
        htmllib.TH("Before", align="CENTER") +
        htmllib.TH("After", align="CENTER") 
        )
        rows.append(x)
        w(htmllib.TABLE("\n".join(rows), border=None, cellpadding=3, cellspacing=0))
        w(htmllib.P())

        w(htmllib.P())
        name = 'Figure 5: This pca plot shows the similarities among your samples after quantile and shiftscale'
        w(htmllib.B(name))
        w(htmllib.P())
       #----------------------------------------------------------------------------
        w(htmllib.HR())
        w(htmllib.A("<methods_normalization>",name="methods_normalization"))
        w(htmllib.CENTER(htmllib.H2("Methods")))
        w(htmllib.H3("1.Quantile File"))
        w('To generate this file, I ran the following analysis:')
        w(htmllib.P())
        for i in range(len(pipelines[0])):
            w('&nbsp&nbsp &nbsp&nbsp &nbsp&nbsp &nbsp&nbsp' +str(i+1)+'. '+pipelines[0][i])
            w(htmllib.P())
        w('I used the following parameters:')
        rows = []
        x = htmllib.TR(
            htmllib.TH("Parameter", align="LEFT") +
            htmllib.TH("Value", align="LEFT") 
            )
        rows.append(x)
        for key in parameters[0].keys():
            x = htmllib.TR(
            htmllib.TD(key, align="LEFT") +
            htmllib.TD(parameters[0][key], align="LEFT") 
            )
            rows.append(x)
        w(htmllib.TABLE("\n".join(rows), border=1, cellpadding=3, cellspacing=0))
        w(htmllib.P())
        #----------------------------------------------------------------------------
        
        w(htmllib.H3("2.Bfrm File"))
        w('To generate this file, I ran the following analysis:')
        w(htmllib.P())
        for i in range(len(pipelines[3])):
            w('&nbsp&nbsp &nbsp&nbsp &nbsp&nbsp &nbsp&nbsp' +str(i+1)+'. '+pipelines[3][i])
            w(htmllib.P())
        w('I used the following parameters:')
        rows = []
        x = htmllib.TR(
            htmllib.TH("Parameter", align="LEFT") +
            htmllib.TH("Value", align="LEFT") 
            )
        rows.append(x)
        for key in parameters[3].keys():
            x = htmllib.TR(
            htmllib.TD(key, align="LEFT") +
            htmllib.TD(parameters[3][key], align="LEFT") 
            )
            rows.append(x)
        w(htmllib.TABLE("\n".join(rows), border=1, cellpadding=3, cellspacing=0))
        w(htmllib.P())
       #----------------------------------------------------------------------------
        w(htmllib.H3("3.Quantile and Combat File"))
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
        #----------------------------------------------------------------------------
        w(htmllib.H3("4.Quantile and Dwd File"))
        w('To generate this file, I ran the following analysis:')
        w(htmllib.P())
        for i in range(len(pipelines[9])):
            w('&nbsp&nbsp &nbsp&nbsp &nbsp&nbsp &nbsp&nbsp' +str(i+1)+'. '+pipelines[9][i])
            w(htmllib.P())
        w('I used the following parameters:')
        rows = []
        x = htmllib.TR(
            htmllib.TH("Parameter", align="LEFT") +
            htmllib.TH("Value", align="LEFT") 
            )
        rows.append(x)
        for key in parameters[9].keys():
            x = htmllib.TR(
            htmllib.TD(key, align="LEFT") +
            htmllib.TD(parameters[9][key], align="LEFT") 
            )
            rows.append(x)
        w(htmllib.TABLE("\n".join(rows), border=1, cellpadding=3, cellspacing=0))
        w(htmllib.P())
        #----------------------------------------------------------------------------
        w(htmllib.H3("5.Quantile and Shiftscale File"))
        w('To generate this file, I ran the following analysis:')
        w(htmllib.P())
        for i in range(len(pipelines[12])):
            w('&nbsp&nbsp &nbsp&nbsp &nbsp&nbsp &nbsp&nbsp' +str(i+1)+'. '+pipelines[12][i])
            w(htmllib.P())
        w('I used the following parameters:')
        rows = []
        x = htmllib.TR(
            htmllib.TH("Parameter", align="LEFT") +
            htmllib.TH("Value", align="LEFT") 
            )
        rows.append(x)
        for key in parameters[12].keys():
            x = htmllib.TR(
            htmllib.TD(key, align="LEFT") +
            htmllib.TD(parameters[12][key], align="LEFT") 
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
        module_utils.write_Betsy_report_parameters_file(
             outfiles,'report.html',starttime,user,jobname)
    except:
        raise 
    finally:
        os.chdir(cwd)
    
    print 'Report:'+ os.path.join(result_folder,'report.html')

