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
        
        data_node1, data_node2 = antecedents
        #write the report.html
        from genomicode import parselib
        from genomicode import htmllib

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
            title = "Clustering Results"
            x = parselib.remove_all_tags(title)
            w(htmllib.HEAD(htmllib.TITLE(x)))
            w("<BODY>")
            w(htmllib.CENTER(htmllib.H1(title)))
            w('I generated a file that contains the expression values of the data set after clustering')
            w(htmllib.P())
            w(htmllib.A(result_files[0], result_files[0]))
            w(htmllib.P())
            w(htmllib.A("Methods", href="#methods_clustering"))
            w(htmllib.P())
            w(htmllib.A(htmllib.IMG(height=500,
                                    src=result_files[1]),
                        href=result_files[1]))
            w(htmllib.P())
            name = 'Figure 1: In this heatmap, each row contains a signature and each column \
            contains a sample from your data set.'

            w(htmllib.B(name))

            w(htmllib.HR())
            w(htmllib.A("<methods_clustering>", name="methods_clustering"))
            w(htmllib.CENTER(htmllib.H2("Methods")))
            w('To generate this file, I ran the following analysis:')
            bie3.plot_network_gv(os.path.join(outfile_folder, "network.png"),
                                 network)
            ##        w(htmllib.P())
            ##        for i in range(len(pipelines[0])):
            ##            w('&nbsp&nbsp &nbsp&nbsp &nbsp&nbsp &nbsp&nbsp' +str(i+1)+'. '+pipelines[0][i])
            ##            w(htmllib.P())
            w(htmllib.A(htmllib.IMG(height=500,
                                    src="network.png"),
                        href="network.png"))
            w(htmllib.P())
            name = 'Figure 1: In this heatmap, each row contains a signature and each column \
            contains a sample from your data set.'

            w('I used the following parameters:')
            w(htmllib.H3("1.Cluster File"))
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
            w(htmllib.H3("2.Cluster heatmap"))
            rows = []
            x = htmllib.TR(htmllib.TH("Parameter",
                                      align="LEFT") + htmllib.TH("Value",
                                                                 align="LEFT"))
            rows.append(x)

            for key in data_node2.data.attributes.keys():
                x = htmllib.TR(htmllib.TD(key,
                                          align="LEFT") +
                               htmllib.TD(data_node2.data.attributes[key],
                                          align="LEFT"))
                rows.append(x)
            w(htmllib.TABLE("\n".join(rows),
                            border=1,
                            cellpadding=3,
                            cellspacing=0))
            w(htmllib.P())
            # Write out the footer.
            #time_str = parselib.pretty_date(time.time())
            #hostname = pretty_hostname()
            w(htmllib.P())
            w(htmllib.HR())
            #w(htmllib.EM(
            #    "This analysis was run on %s on %s. \n" %
            #    (time_str, hostname)))
            w("</BODY>")
            w("</HTML>")
            x = "\n".join(lines) + "\n"
            open(outfile, 'w').write(x)
        except:
            raise



    
    def name_outfile(self, antecedents, user_options):
        filename = 'report'
        return filename


    def find_antecedents(
        self, network, module_id, out_attributes, user_attributes, pool):
        from Betsy import module_utils
        filter1 = module_utils.AntecedentFilter(datatype_name='ClusterFile')
        filter2 = module_utils.AntecedentFilter(datatype_name='Heatmap')
        x = module_utils.find_antecedents(
            network, module_id, user_attributes, pool, filter1, filter2)
        return x
