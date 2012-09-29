#make_diffgenes_report.py
from Betsy import config
from Betsy import module_utils
from Betsy import hash_method
import os
import shutil
from Betsy import protocol_utils
import imghdr
import time
from genomicode import parselib,htmllib
def run(outfiles,parameters,pipelines):
    OUTPUTPATH = config.OUTPUTPATH
    inputid = module_utils.get_inputid(outfiles[0])
    folder_string = hash_method.hash_parameters(
            inputid,pipelines[0],**parameters[0])
    folder_name = 'diffgenes_report_BETSYHASH1_'+folder_string
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
        title = "Differential Expressed Gene Analysis Result"
        x = parselib.remove_all_tags(title)
        w(htmllib.HEAD(htmllib.TITLE(x)))
        w("<BODY>")
        w(htmllib.CENTER(htmllib.H1(title)))
        w('I generated a file that contains the t-test results')
        w(htmllib.P())
        w(htmllib.A(result_files[0], result_files[0]))
        w(htmllib.P())
        w(htmllib.A("Methods",href="#methods_diffgenes"))
        w(htmllib.P())
        w('I generated a file that contains the sam results')
        w(htmllib.P())
        sam_result = os.path.join(result_files[1],'sam_result.txt')
        sam_fig = os.path.join(result_files[1],'sam_plot.png')
        w(htmllib.A(sam_result, sam_result))
        w(htmllib.P())
        #---------------------------------
        name = 'Table 1: Table of significant genes p<0.05 sorted in order of significance'
        w(htmllib.B(name))
        f = file(result_files[0],'rU')
        text = f.readlines()
        f.close()
        header = text[0].split('\t')
        data=[i.split('\t') for i in text[1:]]
        rows = write_table(header,data,50)
        w(htmllib.TABLE("\n".join(rows), border=1, cellpadding=3, cellspacing=0))
        w(htmllib.P())
        if len(text)>51:
            more_genes = len(text)-51
            w(htmllib.A(str(more_genes)+' more genes', result_files[0]))
            w(htmllib.P())
        #-----------------------------------
        w(htmllib.A(htmllib.IMG(height=500,
                src=result_files[2]), href=result_files[2]))
        w(htmllib.P())
        name = 'Figure 1: Heatmap of significant genes'
        w(htmllib.B(name))
        w(htmllib.P())
        #-----------------------------------
        name = 'Table 2: Table of significant annotations'
        w(htmllib.B(name))
        w(htmllib.P())
        f = file(result_files[3],'rU')
        text = f.readlines()
        f.close()
        index = [0,1,2,3,4,5,6,7,9,10]
        header = text[0].split('\t')
        header = [header[i] for i in index]
        data=[i.split('\t') for i in text[1:]]
        data=[[item[i] for i in index] for item in data]
        import math
        p_list=[(data[i][8],i) for i in range(len(data)) if data[i][8]>-math.log(0.05)]
        p_list.sort(reverse=True)
        sort_index = [i[1] for i in p_list]
        data = [data[i] for i in sort_index]
        rows = write_table(header,data,10)
        w(htmllib.TABLE("\n".join(rows), border=1, cellpadding=3, cellspacing=0))
        w(htmllib.P())
        if len(data)>10:
            more_genes = len(data)-10
            w(htmllib.A(str(more_genes)+' more annotations', result_files[3]))
            w(htmllib.P())
        #-----------------------------------
        w(htmllib.A(htmllib.IMG(height=500,
                src=sam_fig), href=sam_fig))
        w(htmllib.P())
        name = 'Figure 2: SAM plot'
        w(htmllib.B(name))
        #-----------------------------------
        w(htmllib.P())
        w('The full result of Gather is in')
        w(htmllib.P())
        w(htmllib.A(result_files[3], result_files[3]))
        w(htmllib.P())
        #-----------------------------------
        w('The result of GSEA is in')
        w(htmllib.P())
        w(htmllib.A(result_files[4], result_files[4]))
        w(htmllib.P())
        w(htmllib.HR())
        #-----------------------------------
        w(htmllib.A("<methods_diffgenes>",name="methods_diffgenes"))
        w(htmllib.CENTER(htmllib.H2("Methods")))
        w(htmllib.H3("1.T-test"))
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
