#make_normalize_report.py
from Betsy import config
from Betsy import module_utils
from Betsy import hash_method
import os
import shutil
from Betsy import protocol_utils
import imghdr
import time
def run(outfiles,parameters,pipelines):
    OUTPUTPATH = config.OUTPUTPATH
    inputid = module_utils.get_inputid(outfiles[0])
    folder_string = hash_method.hash_parameters(
            inputid,pipelines[0],**parameters[0])
    folder_name = 'normalize_report_BETSYHASH1_'+folder_string
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
    module = protocol_utils.import_protocol('normalize_file')
    try:
        lines = []
        w = lines.append
        w("<HTML>")
        title = "Normalization Results"
        x = parselib.remove_all_tags(title)
        w(htmllib.HEAD(htmllib.TITLE(x)))
        w("<BODY>")
        w(htmllib.CENTER(htmllib.H1(title)))
        w('I generated a file that contains the normalized gene expression values')
        w(htmllib.P())
        w(htmllib.A(result_files[0], result_files[0]))
        w(htmllib.P())
        w(htmllib.A("Methods",href="#methods_normalization"))
        w(htmllib.P())
        if pipelines[1] == pipelines[2]:
            w(htmllib.A(htmllib.IMG(height=500,
                src=result_files[1]), href=result_files[1]))
        else:
            rows=[]
            x = htmllib.TR(
                htmllib.TD(htmllib.A(htmllib.IMG(height=500,
                src=result_files[1]), href=result_files[1]), align="CENTER") +
                htmllib.TD(htmllib.A(htmllib.IMG(height=500,
                src=result_files[2]), href=result_files[2]), align="CENTER") 
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
        name = 'Figure 1: This pca plot shows the similarities among your samples'
        w(htmllib.B(name))
        w(htmllib.P())
        w(htmllib.A(htmllib.IMG(height=500,
                src=result_files[3]), href=result_files[3]))
        w(htmllib.P())
        name = 'Figure 2: This boxplot shows the distribution of signal values'
        w(htmllib.B(name))
        w(htmllib.P())
        w(htmllib.A(htmllib.IMG(height=500,
                src=result_files[4]), href=result_files[4]))
        w(htmllib.P())
        name = 'Figure 3: This plot shows the values of ACTB and TUBB genes'
        w(htmllib.B(name))
        w(htmllib.P())
        if len(result_files)==8:
            w(htmllib.A(htmllib.IMG(height=1000,
                src=result_files[5]), href=result_files[5]))
            name = 'Figure 4: This plot shows the value of biotin and housekeeping control genes'
        else:
            w(htmllib.A(htmllib.IMG(height=500,
                src=result_files[5]), href=result_files[5]))
            name = 'Figure 4: This plot shows the average values control genes'
        w(htmllib.P())
        
        
        w(htmllib.B(name))
        if len(result_files)==8:
            w(htmllib.P())
            w(htmllib.A(htmllib.IMG(height=500,
                src=result_files[6]), href=result_files[6]))
            w(htmllib.P())
            name = 'Figure 5: This barplot shows the disctribution control values'
            w(htmllib.B(name))
            w(htmllib.P())
           
        w(htmllib.HR())
        w(htmllib.A("<methods_normalization>",name="methods_normalization"))
        w(htmllib.CENTER(htmllib.H2("Methods")))
        w(htmllib.H3("1.Normalization File"))
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
        w(htmllib.H3("2. PCA analysis"))
        w('I made a principal component plot that shows the similarities among your samples.')
        w(htmllib.P())
        w(htmllib.H3("3. Signal distribution"))
        w('I made a box plot that shows the distribution of signal values.')
        w(htmllib.P())
        w(htmllib.H3("4. Control signal"))
        w('I made two plots that show the values of control signal.')
        w(htmllib.P())
        if len(result_files)==8:
            w(htmllib.H3("5. Control signal"))
            w('I made a bar plot that shows the hybridization controls.')
            w(htmllib.P())
            w('The control file is ')
            w(htmllib.A(result_files[7], result_files[7]))
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
