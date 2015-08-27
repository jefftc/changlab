from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import os
        import shutil
        from Betsy import bie3
        from genomicode import htmllib
        from genomicode import parselib
        outfile_folder = outfile
        outfile = os.path.join(outfile_folder, 'report.html')
        if not os.path.exists(outfile_folder):
            os.mkdir(outfile_folder)
        
        result_files = []
        for data_node in antecedents:
            filename = data_node.identifier
            new_name = os.path.join(outfile_folder, os.path.split(filename)[-1])
            if os.path.isdir(filename):
                shutil.copytree(filename, new_name)
            else:
                shutil.copyfile(filename, new_name)
            result_files.append(os.path.split(new_name)[-1])
        
        data_node1, data_node2, data_node3, data_node4, data_node5 = antecedents

        #write the report.html

        #def highlight(s):
        #    from genomicode import htmllib
        #    return htmllib.SPAN(s, style="background-color:yellow")

        
        #def smaller(s):
        #    from genomicode import htmllib
        #    return htmllib.FONT(s, size=-1)

        
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
            w(htmllib.A("Methods", href="#methods_diffgenes"))
            w(htmllib.P())
            w('I generated a file that contains the sam results')
            w(htmllib.P())
            #sam_result = os.path.join(result_files[1],'sam_result.txt')
            #sam_fig = os.path.join(result_files[1],'sam_plot.png')
            w(htmllib.A(result_files[1], result_files[1]))
            w(htmllib.P())
            #---------------------------------
            name = 'Table 1: Table of significant genes p<0.05 sorted in order of significance'
            w(htmllib.B(name))
            f = file(os.path.join(outfile_folder, result_files[0]), 'rU')
            text = f.readlines()
            f.close()
            header = text[0].split('\t')
            data = [i.split('\t') for i in text[1:]]
            rows = write_table(header, data, 50)
            w(htmllib.TABLE("\n".join(rows),
                            border=1,
                            cellpadding=3,
                            cellspacing=0))
            w(htmllib.P())
            if len(text) > 51:
                more_genes = len(text) - 51
                w(htmllib.A(str(more_genes) + ' more genes', result_files[0]))
                w(htmllib.P())
            #-----------------------------------
            w(htmllib.A(htmllib.IMG(height=500,
                                    src=result_files[2]),
                        href=result_files[2]))
            w(htmllib.P())
            name = 'Figure 1: Heatmap of significant genes'
            w(htmllib.B(name))
            w(htmllib.P())
            #-----------------------------------
            name = 'Table 2: Table of significant annotations'
            w(htmllib.B(name))
            w(htmllib.P())
            f = file(os.path.join(outfile_folder, result_files[3]), 'rU')
            text = f.readlines()
            f.close()
            index = [0, 1, 2, 3, 4, 5, 6, 7, 9, 10]
            header = text[0].split('\t')
            header = [header[i] for i in index]
            data = [i.split('\t') for i in text[1:]]
            data = [[item[i] for i in index] for item in data]
            import math
            p_list = [(data[i][8], i) for i in range(len(data))
                      if data[i][8] > -math.log(0.05)]
            p_list.sort(reverse=True)
            sort_index = [i[1] for i in p_list]
            data = [data[i] for i in sort_index]
            rows = write_table(header, data, 10)
            w(htmllib.TABLE("\n".join(rows),
                            border=1,
                            cellpadding=3,
                            cellspacing=0))
            w(htmllib.P())
            if len(data) > 10:
                more_genes = len(data) - 10
                w(htmllib.A(str(more_genes) + ' more annotations', result_files[3]))
                w(htmllib.P())
            #-----------------------------------
            ##        w(htmllib.A(htmllib.IMG(height=500,
            ##                src=sam_fig), href=sam_fig))
            ##        w(htmllib.P())
            ##        name = 'Figure 2: SAM plot'
            ##        w(htmllib.B(name))
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
            w(htmllib.A("<methods_diffgenes>", name="methods_diffgenes"))
            w(htmllib.CENTER(htmllib.H2("Methods")))
            w(htmllib.H3("1.T-test"))
            w('To generate this file, I ran the following analysis:')
            bie3.plot_network_gv(os.path.join(outfile_folder, "network.png"),
                                 network)
            w(htmllib.P())
            w(htmllib.A(htmllib.IMG(height=500,
                                    src="network.png"),
                        href="network.png"))
            w('I used the following parameters:')
            rows = []
            x = htmllib.TR(htmllib.TH("Parameter",
                                      align="LEFT") + htmllib.TH("Value",
                                                                 align="LEFT"))
            rows.append(x)
            for key in data_node1.data.attributes.keys():
                x = htmllib.TR(htmllib.TD(key,
                                          align="LEFT") +
                               htmllib.TD(data_node1.data.attributes[key],
                                          align="LEFT"))
                rows.append(x)
            w(htmllib.TABLE("\n".join(rows),
                            border=1,
                            cellpadding=3,
                            cellspacing=0))
            w(htmllib.P())
            w("</BODY>")
            w("</HTML>")
            x = "\n".join(lines) + "\n"
            open(outfile, 'w').write(x)
        except:
            raise



    def name_outfile(self, antecedents, user_options):
        filename = 'report'
        return filename


    
    def write_table(header, data, N):
    from genomicode import htmllib
    rows = []
    a = ''
    for i in header:
        a = a + htmllib.TH(i, align='LEFT')
    
    x = htmllib.TR(a)
    rows.append(x)
    for i in range(min(N, len(data))):
        a = ''
        for j in range(len(data[0])):
            a = a + htmllib.TD(data[i][j], align="LEFT")
        x = htmllib.TR(a)
        rows.append(x)
    
    return rows


