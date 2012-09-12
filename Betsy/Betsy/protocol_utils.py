#protocol_utils.py
import os
import shutil
import config
import hash_method
import time
import imghdr
import module_utils

def get_result_folder(protocol,outfiles,parameters,pipeline,foldername):
    OUTPUTPATH = config.OUTPUTPATH
    inputid = module_utils.get_inputid(outfiles[0][0])
    folder_string = hash_method.hash_parameters(
            inputid,pipeline[0][0],**parameters[0][0])
    folder_name = foldername + '_BETSYHASH1_' + folder_string
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
                        
    summarize_report(protocol,final_outfiles,result_folder,parameters,pipeline)
    print 'Report:'+ os.path.join(result_folder,'report.html')
    
def format_prolog_query(
    predicate,parameters,modules):
    parameters=[str(i) for i in parameters]
    str_parameters = ','.join(parameters)
    output = str('[' + str_parameters + '],' + modules)
    query = predicate + '(' + output+')'
    return query

def import_protocol(protocol):
    protocol_name = 'protocols.'+protocol
    mod = __import__(protocol_name)
    mod = getattr(mod, protocol)
    return mod

def pretty_hostname():
    import subprocess
    cmd = "hostname"
    p = subprocess.Popen(
        cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    wh, r = p.stdin, p.stdout
    wh.close()
    hostname = r.read().strip()
    assert hostname, "I could not get the hostname."
    return hostname

def summarize_report(protocol,final_result_files,result_folder,final_parameters,final_pipeline):
    from genomicode import parselib
    from genomicode import htmllib

    def highlight(s):
        return htmllib.SPAN(s, style="background-color:yellow")
    def smaller(s):
        return htmllib.FONT(s, size=-1)
    cwd = os.getcwd()
    os.chdir(result_folder)
    module = import_protocol(protocol)  
    try:
        lines = []
        w = lines.append
        w("<HTML>")
        title = "Pipeline Report"
        
        x = parselib.remove_all_tags(title)
        w(htmllib.HEAD(htmllib.TITLE(x)))
        w("<BODY>")
        w(htmllib.CENTER(htmllib.H1(title)))
        for i in range(len(final_result_files)):
            for j in range(len(final_result_files[i])):
                result_files = final_result_files[i][j]
                parameters = final_parameters[i][j]
                pipeline = final_pipeline[i][j]
                #w(htmllib.H3("I.  Parameters"))

                # Make a table with each parameter.
                rows = []
                x = htmllib.TR(
                    htmllib.TH("Parameter", align="LEFT") +
                    htmllib.TH("Value", align="LEFT") 
                    )
                rows.append(x)
                for key in parameters.keys():
                    x = htmllib.TR(
                    htmllib.TD(key, align="LEFT") +
                    htmllib.TD(parameters[key], align="LEFT") 
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
                w('--->'.join(pipeline))
                w(htmllib.P())
                #w(htmllib.H3("III.  Results"))
                output_type = module.OUTPUTS
                result = result_files
                flag = False
                j = (j+1)%len(module.report)-1
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
                            w(module.report[j])
                            w(htmllib.P())
                        else:
                            flag = True
                    else:
                        flag = True
                if flag:
                    
                    w(module.report[j]+'is shown in: %s'
                        % htmllib.A(prob_file, href=prob_file))
                   
        # Write out the footer.
        time_str = parselib.pretty_date(time.time())
        hostname = pretty_hostname()
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
        
def check_parameters(parameters):
    allow_type = ['string','float','integer','list']
    for value in parameters.itervalues():
        if not isinstance(value,list):
            assert value in allow_type,'%s is not a allow type'%value
    return True

def check_default(default,parameters):
    for key in default.keys():
        assert key in parameters.keys(),'%s is not a valid parameter'%key
        value = str(default[key])
        if isinstance(parameters[key],list):
            assert value in parameters[key],'%s is not correct value for %s'%(value,key)
        elif parameters[key] == 'float':
            assert module_utils.isnumber(value),'%s is not a float number for %s'%(value,key)
        elif parameters[key] == 'string':
            assert isinstance(value,str),'%s is not a string for %s'%(str(value),key)
        elif parameters[key] == 'integer':
            assert value.isdigit(),'%s is not a digit for %s'%(str(value),key)
        elif parameters[key] == 'list':
            assert isinstance(value,str),'%s is not a string which will\
                                convert to a list for %s'%(str(value),key)
    return True

