#protocol_utils.py
import os
import shutil
import Betsy_config
import hash_method
import time
import imghdr
import module_utils

##def get_result_folder(protocol,outfiles,parameters,pipeline):
##    OUTPUTPATH = Betsy_config.OUTPUTPATH
##    filename = os.path.split(outfiles[0][0])[-1]
##    if '_BETSYHASH1_' in filename: 
##        inputid = '_'.join(filename.split('_')[:-2])
##    else:
##        inputid = filename
##    for i in range(len(outfiles[0])):
##        result_files = []
##        folder_string = hash_method.hash_parameters(
##            inputid,pipeline[0][i],**parameters[0][i])
##        folder_name = 'result_folder_BETSYHASH1_'+folder_string
##        result_folder = os.path.join(OUTPUTPATH,folder_name)
##        if not os.path.exists(result_folder):
##            os.mkdir(result_folder)
##        for j in range(len(outfiles)):
##            if len(outfiles[j]) == 1:
##                final_output = os.path.split(outfiles[j][0])[-1]
##                result_file = os.path.join(result_folder,final_output)
##                result_files.append(outfiles[j][0])
##                if not os.path.exists(result_file):
##                    if os.path.isdir(outfiles[j][0]):
##                        shutil.copytree(outfiles[j][0],result_file)
##                    else:
##                        shutil.copyfile(outfiles[j][0],result_file) 
##            elif len(outfiles[j]) > 1:
##                final_output = os.path.split(outfiles[j][i])[-1]
##                result_file = os.path.join(result_folder,final_output)
##                result_files.append(outfiles[j][i])
##                if not os.path.exists(result_file):
##                    if os.path.isdir(outfiles[j][i]):
##                        shutil.copytree(outfiles[j][i],result_file)
##                    else:
##                        shutil.copyfile(outfiles[j][i],result_file)
##            elif len(outfiles[j]) == 0:
##                result_files.append(None)
##        summarize_report(protocol,result_files,result_folder,parameters[0][i],pipeline[0][i])
def get_result_folder(protocol,outfiles,parameters,pipeline):
    OUTPUTPATH = Betsy_config.OUTPUTPATH
    filename = os.path.split(outfiles[0][0])[-1]
    if '_BETSYHASH1_' in filename: 
        inputid = '_'.join(filename.split('_')[:-2])
    else:
        inputid = filename
    for i in range(len(outfiles[0])):
        result_files = []
        folder_string = hash_method.hash_parameters(
            inputid,pipeline[0][i],**parameters[0][i])
        folder_name = 'result_folder_BETSYHASH1_'+folder_string
        result_folder = os.path.join(OUTPUTPATH,folder_name)
        if not os.path.exists(result_folder):
            os.mkdir(result_folder)
        for j in range(len(outfiles)):
            if len(outfiles[j]) == 1:
                final_output = os.path.split(outfiles[j][0])[-1]
                result_file = os.path.join(result_folder,final_output)
                result_files.append(outfiles[j][0])
                if not os.path.exists(result_file):
                    if os.path.isdir(outfiles[j][0]):
                        shutil.copytree(outfiles[j][0],result_file)
                    else:
                        shutil.copyfile(outfiles[j][0],result_file) 
            elif len(outfiles[j]) > 1:
                final_output = os.path.split(outfiles[j][i])[-1]
                result_file = os.path.join(result_folder,final_output)
                result_files.append(outfiles[j][i])
                if not os.path.exists(result_file):
                    if os.path.isdir(outfiles[j][i]):
                        shutil.copytree(outfiles[j][i],result_file)
                    else:
                        shutil.copyfile(outfiles[j][i],result_file)
            elif len(outfiles[j]) == 0:
                result_files.append(None)
        summarize_report(protocol,result_files,result_folder,parameters[0][i],pipeline[0][i])
def format_prolog_query(
    predicate,parameters,modules):
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

def summarize_report(protocol,result_files,result_folder,parameters,pipeline):
    import subprocess
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
        title = "Pipeline Report"
        
        x = parselib.remove_all_tags(title)
        w(htmllib.HEAD(htmllib.TITLE(x)))
        w("<BODY>")
        w(htmllib.CENTER(htmllib.H1(title)))

        w(htmllib.H3("I.  Parameters"))

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
        w(htmllib.H3("II.  Pipeline Sequence"))
        w(htmllib.P())
        w("The modules run in the follow sequence:")
        w(htmllib.P())
        w(htmllib.P())
        w('--->'.join(pipeline))
        w(htmllib.P())
        #show all the figures
        module = import_protocol(protocol)
        all_description = {
            'cluster_heatmap':"In this heatmap, each row contains a signature and each column "
                               "contains a sample from your data set.",
             'cluster_file' : 'the expression value of the data set after clustering',
            'signal_file' : 'the expression value of the data set after normalization',
            'pca_plot_in' : 'the pca plot of the data set before normalization',
            'pca_plot_out' : 'the pca plot of the data set after normalization',
            'intensity_plot': 'the intersity of the signal in the data set after normalization',
            'biotin_plot': 'the value of biotin and housekeeping in different sample in the control file',
            'actb_plot': 'the value of ACTB and TUBB in different sample before normalization in the data set',
            'hyb_bar_plot': 'the value of hybridization controls in the data set',
            'svm_predictions':'the svm predictions',
            'weightedVoting':'the results file of weighted Voting ',
            'loocv': ' the result of the leave one out cross validation',
            'differential_expressed_genes':'the result of the differential_expressed_genes ',
            'class_neighbors':'the result of the class_neighbors'
            }
        #put the image in the same folder
        w(htmllib.H3("III.  Results"))
        output_type = module.OUTPUTS
        for i in range(len(output_type)):
            result = result_files[i]
            flag = False
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
                        w(all_description[output_type[i]])
                        w(htmllib.P())
                    else:
                        flag = True
                else:
                    flag = True
            if flag:
                w(all_description[output_type[i]]+'is shown in: %s'
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

