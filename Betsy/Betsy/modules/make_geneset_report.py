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
        
        data_node1, data_node2 = antecedents
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
            title = "Geneset Analysis Results"
            x = parselib.remove_all_tags(title)
            w(htmllib.HEAD(htmllib.TITLE(x)))
            w("<BODY>")
            w(htmllib.CENTER(htmllib.H1(title)))
            w('I generated a file that contains the analysis result of the geneset')
            w(htmllib.P())
            w(htmllib.A(result_files[0], result_files[0]))
            w(htmllib.P())
            w(htmllib.A("Methods", href="#methods"))
            w(htmllib.P())
            filenames = os.listdir(os.path.join(outfile_folder, result_files[1]))
            c = 0
            for filename in filenames:
                c = c + 1
                w(htmllib.A(
                    htmllib.IMG(height=500,
                                src=os.path.join(result_files[1], filename)),
                    href=os.path.join(result_files[1], filename)))
                w(htmllib.P())
                name = 'Figure ' + str(c) + ': Geneset Plot.'
                w(htmllib.B(name))
            w(htmllib.HR())
            w(htmllib.A("<methods>", name="methods"))
            w(htmllib.CENTER(htmllib.H2("Methods")))
            w(htmllib.H3("1.Result File"))
            w('To generate this file, I ran the following analysis:')
            bie3.plot_network_gv(os.path.join(outfile_folder, "network.png"),
                                 network)
            w(htmllib.A(htmllib.IMG(height=500,
                                    src="network.png"),
                        href="network.png"))
            w(htmllib.P())

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
