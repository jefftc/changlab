#make_geneset_report.py
from Betsy import config
from Betsy import module_utils
from Betsy import hash_method
import os
import shutil
import imghdr
import time
from genomicode import parselib
from genomicode import htmllib


def run(outfiles, parameters, pipelines):
    OUTPUTPATH = config.OUTPUTPATH
    inputid = module_utils.get_inputid(outfiles[0])
    folder_string = hash_method.hash_parameters(
            inputid, pipelines[0], **parameters[0])
    folder_name = 'geneset_report_BETSYHASH1_' + folder_string
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
