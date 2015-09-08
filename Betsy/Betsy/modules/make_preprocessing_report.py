from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import os
        import shutil
        #import time
        from genomicode import parselib
        from genomicode import htmllib
        from Betsy import bie3
        
        out_path = outfile
        outfile = os.path.join(out_path, 'report.html')
        if not os.path.exists(out_path):
            os.mkdir(out_path)

        (data_node1, data_node2, data_node3, data_node4,
         data_node5, data_node6) = antecedents
        
        # Make a list of the (relative) files for each input nodes:
        # 0  SignalFile     Preprocessed gene expression data.
        # 1  IntensityPlot  Box plot of signal intensity values.
        # 2  ControlPlot    AFFX control probes.
        # 3  PcaPlot        Has normalization.
        # 4  ActbPlot       From _SignalFile_Impute.
        # 5  PcaPlot        No normalization or anything.
        result_files = []
        for index, data_node in enumerate(antecedents):
            in_filename = data_node.identifier   # full path
            in_path, in_file = os.path.split(in_filename)
            out_file = in_file
            #rename one of the pcaplot filename
            if index == 3:
                out_file = "after_%s" % in_file
            out_filename = os.path.join(out_path, out_file)
            if os.path.isdir(in_filename):
                shutil.copytree(in_filename, out_filename)
            else:
                shutil.copyfile(in_filename, out_filename)
            result_files.append(out_file)
        
        #write the report.html

        IMG_HEIGHT = 400

        lines = []
        w = lines.append
        
        w("<HTML>")
        title = "Normalization Results"
        x = parselib.remove_all_tags(title)
        w(htmllib.HEAD(htmllib.TITLE(x)))
        w("<BODY>")
        w(htmllib.CENTER(htmllib.H1(title)))

        # Provide a link to the signal file.
        w('Preprocessed signal values: ')
        # TODO: Show the attributes for this data.
        w(htmllib.A(result_files[0], result_files[0]))
        w(htmllib.P())
        
        ##w(htmllib.A("Methods", href="#methods_normalization"))
        ##w(htmllib.P())
        ##        if pipelines[1] == pipelines[2]:
        ##            w(htmllib.A(htmllib.IMG(height=500,
        ##                src=result_files[1]), href=result_files[1]))
        ##        else:

        # Show the PCA plot before and after normalization.
        rows = []
        col1 = htmllib.A(
            htmllib.IMG(height=IMG_HEIGHT, src=result_files[5]),
            href=result_files[5])
        col2 = htmllib.A(
            htmllib.IMG(height=IMG_HEIGHT, src=result_files[3]),
            href=result_files[3])
        x = htmllib.TR(htmllib.TD(col1, align="CENTER") +
                       htmllib.TD(col2, align="CENTER"))
        rows.append(x)
        x = htmllib.TR(htmllib.TH("Before", align="CENTER") +
                       htmllib.TH("After", align="CENTER"))
        rows.append(x)
        w(htmllib.TABLE(
            "\n".join(rows), border=None, cellpadding=3, cellspacing=0))
        w(htmllib.P())
        w(htmllib.P())
        name = 'Figure 1: This pca plot shows the similarities among your samples'
        w(htmllib.B(name))
        w(htmllib.P())

        # Show the distribution of the signal values.
        w(htmllib.A(htmllib.IMG(
            height=IMG_HEIGHT, src=result_files[1]), href=result_files[1]))
        w(htmllib.P())
        name = 'Figure 2: This boxplot shows the distribution of signal values'
        w(htmllib.B(name))
        w(htmllib.P())

        # Show the actin and tubulin values.
        w(htmllib.A(htmllib.IMG(
            height=IMG_HEIGHT, src=result_files[4]), href=result_files[4]))
        w(htmllib.P())
        name = 'Figure 3: This plot shows the values of ACTB and TUBB genes'
        w(htmllib.B(name))
        w(htmllib.P())

        # Affymetrix control genes.
        w(htmllib.A(htmllib.IMG(
            height=IMG_HEIGHT, src=result_files[2]), href=result_files[2]))
        name = 'Figure 4: This plot shows the average values Affymetrix control genes'
        w(htmllib.P())
        w(htmllib.B(name))

        w(htmllib.HR())
        w(htmllib.A("<methods_normalization>", name="methods_normalization"))
        w(htmllib.CENTER(htmllib.H2("Methods")))
        w(htmllib.H3("1.Normalization File"))
        w('To generate this file, I ran the following analysis:')
        w(htmllib.P())
        bie3.plot_network_gv(
            os.path.join(out_path, "network.png"), network)
        w(htmllib.A(htmllib.IMG(
            height=IMG_HEIGHT, src="network.png"), href="network.png"))
        w(htmllib.P())

        w('I used the following parameters:')
        rows = []
        x = htmllib.TR(htmllib.TH(
            "Parameter", align="LEFT") + htmllib.TH("Value", align="LEFT"))
        rows.append(x)
        for key in data_node1.data.attributes.keys():
            x = htmllib.TR(htmllib.TD(
                key, align="LEFT") + htmllib.TD(
                               data_node1.data.attributes[key],
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


    def name_outfile(self, antecedents, user_options):
        return 'report'
