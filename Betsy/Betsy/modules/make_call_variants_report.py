#make_call_variants_report.py

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
    folder_name = 'call_variants_report_BETSYHASH1_'+folder_string
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
        title = "Calling Variants Results"
        x = parselib.remove_all_tags(title)
        w(htmllib.HEAD(htmllib.TITLE(x)))
        w("<BODY>")
        w(htmllib.CENTER(htmllib.H1(title)))
        w('I generated a file that contains the variants by GATK standard method')
        w(htmllib.P())
        w(htmllib.A(result_files[0], result_files[0]))
        w(htmllib.P())
        w(htmllib.A("Methods",href="#methods_calling_variants"))
        w(htmllib.P())
        #------------------------------------------------------------------------
        w('I generated a file that contains the variants by mpileup method')
        w(htmllib.P())
        w(htmllib.A(result_files[1], result_files[1]))
        w(htmllib.P())
        #------------------------------------------------------------------------
        w(htmllib.HR())
        w(htmllib.A("<methods_calling_variants>",name="methods_calling_variants"))
        w(htmllib.CENTER(htmllib.H2("Methods")))
        w(htmllib.H3("1.GATK Calling Variants"))
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
        w(htmllib.H3("2.Mpileup Calling Variants"))
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
        #----------------------------------------------------------------------------
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

