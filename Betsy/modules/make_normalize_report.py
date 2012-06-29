#make_normalize_report.py
import Betsy_config
import module_utils
import hash_method
import os
import shutil
import protocol_utils
import imghdr
import time
def run(outfiles,parameters,pipeline):
    OUTPUTPATH = Betsy_config.OUTPUTPATH
    inputid = module_utils.get_inputid(outfiles[0][0])
    folder_string = hash_method.hash_parameters(
            inputid,pipeline[0][0],**parameters[0][0])
    folder_name = 'normalize_report_BETSYHASH1_'+folder_string
    result_folder = os.path.join(OUTPUTPATH,folder_name)
    if not os.path.exists(result_folder):
        os.mkdir(result_folder)
    final_outfiles = []
    for i in range(len(outfiles)):
        result_files = []
        new_outfiles = []
        for j in range(len(outfiles[i])):
                final_output = os.path.split(outfiles[i][j])[-1]
                output_folder = os.path.split(os.path.split(outfiles[i][j])[-2])[-1]
                new_filename = output_folder+'_'+ final_output
                new_outfiles.append(new_filename)
                result_file = os.path.join(result_folder,new_filename)
                if not os.path.exists(result_file):
                    if os.path.isdir(outfiles[i][j]):
                        shutil.copytree(outfiles[i][j],result_file)
                    else:
                        shutil.copyfile(outfiles[i][j],result_file)
        final_outfiles.append(new_outfiles)
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
        title = "Pipeline Report"
        
        x = parselib.remove_all_tags(title)
        w(htmllib.HEAD(htmllib.TITLE(x)))
        w("<BODY>")
        w(htmllib.CENTER(htmllib.H1(title)))
        for i in range(len(final_outfiles)):
            for j in range(len(final_outfiles[i])):
                result_files = final_outfiles[i][j]
                parameter = parameters[i][j]
                single_pipeline = pipeline[i][j]
                #w(htmllib.H3("I.  Parameters"))
                # Make a table with each parameter.
                rows = []
                x = htmllib.TR(
                    htmllib.TH("Parameter", align="LEFT") +
                    htmllib.TH("Value", align="LEFT") 
                    )
                rows.append(x)
                for key in parameter.keys():
                    x = htmllib.TR(
                    htmllib.TD(key, align="LEFT") +
                    htmllib.TD(parameter[key], align="LEFT") 
                    )
                    rows.append(x)
                w(htmllib.TABLE("\n".join(rows), border=1, cellpadding=3, cellspacing=0))
                w(htmllib.P())
                w(htmllib.B("Table 1: Parameters set up."))
                w("The parameters of the result files are shown above")
                w(htmllib.P())
               #show the pipeline sequence
                #w(htmllib.H3("II.  Pipeline Sequence"))
                w(htmllib.P())
                w("The modules run in the follow sequence:")
                w(htmllib.P())
                w(htmllib.P())
                w('--->'.join(single_pipeline))
                w(htmllib.P())
                #w(htmllib.H3("III.  Results"))
                output_type = module.OUTPUTS
                result = result_files
                flag = False
                if parameters[0][0]['preprocess'] == 'illumina':
                    report = module.report
                else:
                    report = module.report[0:4]
                    
                j = (j+1)%len(report)-1
                    
                if result:
                    prob_file = os.path.split(result)[-1]
                    if not(os.path.isdir(prob_file)):
                        if imghdr.what(prob_file) == 'png':
                            w(htmllib.P())
                            w(htmllib.A(htmllib.IMG(height=500,
                                    src=prob_file), href=prob_file))
                            w(htmllib.P())
                            name = 'Figure:' + output_type[i]
                            w(htmllib.B(name))
                            w(report[j])
                            w(htmllib.P())
                        else:
                            flag = True
                    else:
                        flag = True
                if flag:
                    
                    w(report[j]+'is shown in: %s'
                        % htmllib.A(prob_file, href=prob_file))
        # Write out the footer.
        time_str = parselib.pretty_date(time.time())
        hostname = protocol_utils.pretty_hostname()
        w(htmllib.P())
        w(htmllib.HR())
        w(htmllib.EM(
            "This analysis was run on %s on %s. \n" %
            (time_str, hostname)))
        w("</BODY>")
        w("</HTML>")
        x = "\n".join(lines) + "\n"
        open('report.html', 'w').write(x)
    except:
        raise 
    finally:
        os.chdir(cwd)              
    
    print 'Report:'+ os.path.join(result_folder,'report.html')
