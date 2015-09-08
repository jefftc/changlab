from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        outfile):
        import os
        from genomicode import parselib
        from genomicode import htmllib
        from Betsy import reportlib
        from Betsy import bie3
        
        out_path = outfile
        outfile = os.path.join(out_path, 'report.html')
        if not os.path.exists(out_path):
            os.mkdir(out_path)
        
        #(data_node1, data_node2, data_node3, data_node4, data_node5,
        # data_node6, data_node7, data_node8, data_node9) = antecedents

        # Make a list of the (relative) files for each input nodes:
        # 0  SignalFile         Preprocessed gene expression data.
        # 1  IntensityPlot      Box plot of signal intensity values.
        # 2  BiotinPlot
        # 3  PcaPlot            Has normalization.
        # 4  ActbPlot           From _SignalFile_Impute.
        # 5  PcaPlot            No normalization or anything.
        # 6  HousekeepingPlot
        # 7  Hyb_barPlot        No normalization or anything.
        # 8  ControlFile        File with Illumina control probes.

        def rename_pca_file(i, in_file):
            out_file = in_file
            if i == 3:
                out_file = 'after_%s' % in_file
            return out_file
                
        filenames = reportlib.extract_filenames(
            antecedents, out_path, rename_pca_file)
        for x in filenames:
            in_file, out_file, in_filename, out_filename = x
            reportlib.copy_file_or_path(in_filename, out_filename)
        signal_file = filenames[0][1]
        intensity_file = filenames[1][1]
        biotin_file = filenames[2][1]
        after_pca_file = filenames[3][1]
        actb_file = filenames[4][1]
        before_pca_file = filenames[5][1]
        housekeeping_file = filenames[6][1]
        hyb_file = filenames[7][1]
        control_file = filenames[8][1]

        signal_node = antecedents[0]
        
        
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
        w(htmllib.A(signal_file, signal_file))
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
            htmllib.IMG(height=IMG_HEIGHT, src=before_pca_file),
            href=before_pca_file)
        col2 = htmllib.A(
            htmllib.IMG(height=IMG_HEIGHT, src=after_pca_file),
            href=after_pca_file)
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
        w(htmllib.A(htmllib.IMG(height=IMG_HEIGHT, src=intensity_file),
                    href=intensity_file))
        w(htmllib.P())
        name = 'Figure 2: This boxplot shows the distribution of signal values'
        w(htmllib.B(name))
        w(htmllib.P())

        
        w(htmllib.A(htmllib.IMG(height=IMG_HEIGHT, src=actb_file),
                    href=actb_file))
        w(htmllib.P())
        name = 'Figure 3: This plot shows the values of ACTB and TUBB genes'
        w(htmllib.B(name))
        w(htmllib.P())

        w(htmllib.A(htmllib.IMG(height=IMG_HEIGHT, src=biotin_file),
                    href=biotin_file))
        w(htmllib.P())
        x = 'Figure 4: This plot shows the value of biotin control genes'
        w(htmllib.B(x))
        w(htmllib.P())

        w(htmllib.A(htmllib.IMG(height=IMG_HEIGHT, src=housekeeping_file),
                    href=housekeeping_file))
        w(htmllib.P())
        x = 'Figure 5: This plot shows the value of housekeeping control genes'
        w(htmllib.B(x))
        w(htmllib.P())

        
        w(htmllib.A(htmllib.IMG(height=IMG_HEIGHT, src=hyb_file),
                    href=hyb_file))
        w(htmllib.P())
        x = 'Figure 6: This barplot shows the distribution control values'
        w(htmllib.B(x))
        w(htmllib.P())


        # Methods.
        w(htmllib.A("<methods_normalization>", name="methods_normalization"))
        w(htmllib.CENTER(htmllib.H2("Methods")))

        w(htmllib.H3("1.Normalization File"))
        w('To generate this file, I ran the following analysis:')
        w(htmllib.P())
        bie3.plot_network_gv(os.path.join(out_path, "network.png"),
                             network)
        w(htmllib.A(
            htmllib.IMG(height=IMG_HEIGHT, src="network.png"),
            href="network.png"))
        w(htmllib.P())


        w('I used the following parameters:')
        rows = []
        x = htmllib.TR(htmllib.TH("Parameter",
                                  align="LEFT") + htmllib.TH("Value",
                                                             align="LEFT"))
        rows.append(x)
        for key in signal_node.data.attributes.keys():
            x = htmllib.TR(htmllib.TD(key,
                                      align="LEFT") +
                           htmllib.TD(signal_node.data.attributes[key],
                                      align="LEFT"))
            rows.append(x)
        w(htmllib.TABLE("\n".join(rows),
                        border=1,
                        cellpadding=3,
                        cellspacing=0))
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
        w(htmllib.H3("5. Control signal"))
        w('I made a bar plot that shows the hybridization controls.')
        w(htmllib.P())
        w('The control file is ')
        w(htmllib.A(control_file, control_file))
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
        filename = 'report'
        return filename


    
