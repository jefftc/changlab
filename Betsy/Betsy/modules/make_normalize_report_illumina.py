#make_normalize_report_illumina.py

import os
import shutil
import time
from Betsy import bie3
from Betsy import rulebase
from Betsy import config
from Betsy import module_utils
from Betsy import hash_method


def run(network, antecedents, out_attributes, user_options, num_cores):
    outfile_folder = name_outfile(antecedents, user_options)
    outfile = os.path.join(outfile_folder, 'report.html')
    if not os.path.exists(outfile_folder):
        os.mkdir(outfile_folder)
    result_files = []
    for index, data_node in enumerate(antecedents):
        filename = data_node.identifier
        new_name = os.path.join(outfile_folder, os.path.split(filename)[-1])
        if index == 2:
            new_name = os.path.join(outfile_folder,
                                    'after_' + os.path.split(filename)[-1])
        if os.path.isdir(filename):
            shutil.copytree(filename, new_name)
        else:
            shutil.copyfile(filename, new_name)
        result_files.append(os.path.split(new_name)[-1])
    (data_node1, data_node2, data_node3, data_node4, data_node5, data_node6,
     data_node7, data_node8, data_node9) = antecedents
    #write the report.html
    from genomicode import parselib
    from genomicode import htmllib

    def highlight(s):
        return htmllib.SPAN(s, style="background-color:yellow")

    def smaller(s):
        return htmllib.FONT(s, size=-1)

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
        w(htmllib.A("Methods", href="#methods_normalization"))
        w(htmllib.P())
        ##        if pipelines[1] == pipelines[2]:
        ##            w(htmllib.A(htmllib.IMG(height=500,
        ##                src=result_files[1]), href=result_files[1]))
        ##        else:
        rows = []
        x = htmllib.TR(htmllib.TD(htmllib.A(htmllib.IMG(height=500,
                                                        src=result_files[1]),
                                            href=result_files[1]),
                                  align="CENTER") +
                       htmllib.TD(htmllib.A(htmllib.IMG(height=500,
                                                        src=result_files[2]),
                                            href=result_files[2]),
                                  align="CENTER"))
        rows.append(x)
        x = htmllib.TR(htmllib.TH("Before",
                                  align="CENTER") + htmllib.TH("After",
                                                               align="CENTER"))
        rows.append(x)
        w(htmllib.TABLE("\n".join(rows),
                        border=None,
                        cellpadding=3,
                        cellspacing=0))
        w(htmllib.P())

        w(htmllib.P())
        name = 'Figure 1: This pca plot shows the similarities among your samples'
        w(htmllib.B(name))
        w(htmllib.P())
        w(htmllib.A(htmllib.IMG(height=500,
                                src=result_files[3]),
                    href=result_files[3]))
        w(htmllib.P())
        name = 'Figure 2: This boxplot shows the distribution of signal values'
        w(htmllib.B(name))
        w(htmllib.P())
        w(htmllib.A(htmllib.IMG(height=500,
                                src=result_files[4]),
                    href=result_files[4]))
        w(htmllib.P())
        name = 'Figure 3: This plot shows the values of ACTB and TUBB genes'
        w(htmllib.B(name))
        w(htmllib.P())

        w(htmllib.A(htmllib.IMG(height=500,
                                src=result_files[5]),
                    href=result_files[5]))
        name = 'Figure 4: This plot shows the value of biotin control genes'
        w(htmllib.P())
        w(htmllib.B(name))
        w(htmllib.P())
        w(htmllib.A(htmllib.IMG(height=500,
                                src=result_files[6]),
                    href=result_files[6]))
        name = 'Figure 5: This plot shows the value of  housekeeping control genes'
        w(htmllib.HR())
        w(htmllib.P())
        w(htmllib.B(name))
        w(htmllib.P())
        w(htmllib.A(htmllib.IMG(height=500,
                                src=result_files[7]),
                    href=result_files[7]))
        w(htmllib.P())
        name = 'Figure 6: This barplot shows the disctribution control values'
        w(htmllib.B(name))
        w(htmllib.P())

        w(htmllib.A("<methods_normalization>", name="methods_normalization"))
        w(htmllib.CENTER(htmllib.H2("Methods")))

        w(htmllib.H3("1.Normalization File"))
        w('To generate this file, I ran the following analysis:')
        w(htmllib.P())
        bie3.plot_network_gv(os.path.join(outfile_folder, "network.png"),
                             network)
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
        w(htmllib.A(result_files[8], result_files[8]))
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
        open(outfile, 'w').write(x)
    except:
        raise
    out_node = bie3.Data(rulebase.NormalizeReportFile, **out_attributes)
    out_object = module_utils.DataObject(out_node, outfile)
    return out_object


def name_outfile(antecedents, user_options):
    filename = 'report'
    outfile = os.path.join(os.getcwd(), filename)
    return outfile


def get_out_attributes(antecedents, out_attributes):
    return out_attributes


def make_unique_hash(pipeline, antecedents, out_attributes, user_options):
    (data_node1, data_node2, data_node3, data_node4, data_node5, data_node6,
     data_node7, data_node8, data_node9) = antecedents
    identifier = data_node1.identifier
    return module_utils.make_unique_hash(identifier, pipeline, out_attributes,
                                         user_options)


def find_antecedents(network, module_id, out_attributes, user_attributes,
                     pool):
    data_node1 = module_utils.get_identifier(network, module_id, pool,
                                             user_attributes,
                                             datatype='SignalFile')
    data_node2, data_node3 = module_utils.find_pcaplots(network, pool,
                                                        module_id)
    data_node4 = module_utils.get_identifier(network, module_id, pool,
                                             user_attributes,
                                             datatype='IntensityPlot')
    data_node5 = module_utils.get_identifier(network, module_id, pool,
                                             user_attributes,
                                             datatype='ActbPlot')
    data_node6 = module_utils.get_identifier(network, module_id, pool,
                                             user_attributes,
                                             datatype='BiotinPlot')
    data_node7 = module_utils.get_identifier(network, module_id, pool,
                                             user_attributes,
                                             datatype='HousekeepingPlot')
    data_node8 = module_utils.get_identifier(network, module_id, pool,
                                             user_attributes,
                                             datatype='Hyb_barPlot')
    data_node9 = module_utils.get_identifier(network, module_id, pool,
                                             user_attributes,
                                             datatype='ControlFile')
    return (data_node1, data_node2, data_node3, data_node4, data_node5,
            data_node6, data_node7, data_node8, data_node9)
