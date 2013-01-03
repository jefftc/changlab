#!/usr/bin/env python
import os
import argparse
import numpy
import sys
import arrayio
from genomicode import filelib, jmath, genesetlib, config
from genomicode import matrixlib, Matrix, colorlib


def guess_file_type(filename):
    M = arrayio.read(filename)
    headers = M._row_order
    if headers == ['FILE', 'SAMPLE']:
        return 'geneset_file'
    else:
        return 'gene_expression_file'


def calc_km(survival, dead, group):
    R = jmath.start_R()
    jmath.R_equals(survival, 'survival')
    jmath.R_equals(dead, 'dead')
    new_group = ['"' + i + '"' for i in group]
    jmath.R_equals(new_group, 'group')
    R('x <- calc.km.multi(survival, dead, group)')
    c = R['x']
    p_value = c.rx2('p.value')[0]
    unique_group = list(set(group))
    unique_group.sort()
    surv90 = [''] * len(unique_group)
    surv50 = [''] * len(unique_group)
    for k in range(len(unique_group)):
        surv90[k] = c.rx2('surv').rx2(unique_group[k]).rx2('surv.90')[0]
        surv50[k] = c.rx2('surv').rx2(unique_group[k]).rx2('surv.50')[0]
    MAX_SURV = 1e10
    med_high, med_low = surv50[-1], surv50[0]
    direction = ''
    if (str(med_high), str(med_low)) == ('NA', 'NA'):
        med_high, med_low = surv90[-1], surv90[0]
    if str(med_high) == 'NA':
        med_high = MAX_SURV
    if str(med_low) == 'NA':
        med_low = MAX_SURV
    assert med_low <= MAX_SURV and med_high <= MAX_SURV
    if med_high > med_low:
        direction = ('Low expression has shorter survival.')
    elif med_high < med_low:
        direction = ('High expression has shorter survival.')
    if p_value >= 0.05:
        direction = ''
    surv90 = ['' if (str(k) == 'NA') else str(k) for k in surv90]
    surv50 = ['' if (str(k) == 'NA') else str(k) for k in surv50]
    return p_value, surv90, surv50, direction


def parse_cutoffs(cutoffs):
    if cutoffs:
        cutoffs = cutoffs.split(',')
        cutoffs = [float(i) for i in cutoffs]
        cutoffs.sort()
    if not cutoffs:
        cutoffs = [0.5]
    for cutoff in cutoffs:
        assert cutoff > 0 and cutoff < 1.0, 'cutoff should be between 0 and 1'
    return cutoffs


def parse_genes(genes):
    if genes:
        gene_list = []
        for i in genes:
            gene_list.extend(i.split(','))
        return gene_list
    else:
        return None


def intersect_clinical_matrix_sample_name(M, clinical_data):
    sample_name = M._col_names['_SAMPLE_NAME']
    name_in_clinical = None
    clinical_dict = dict()
    for g in clinical_data:
        name, description, values = g
        clinical_dict[name] = values
    for g in clinical_dict:
        n = len(list(set(sample_name).intersection(set(clinical_dict[g]))))
        if n > 0:
            name_in_clinical = clinical_dict[g]
    assert name_in_clinical, ('cannot match the sample name in '
                              'clinical data and gene data')
    return clinical_dict, name_in_clinical


def convert_genesetfile2matrix(filename):
    matrix = [x for x in filelib.read_cols(filename)]
    matrix = jmath.transpose(matrix)
    col_names = {}
    col_names['_SAMPLE_NAME'] = matrix[1][1:]
    row_names = {}
    row_names['geneset'] = []
    data = []
    for line in matrix[2:]:
        single_data = [float(i) for i in line[1:]]
        data.append(single_data)
        row_names['geneset'].append(line[0])
    M = Matrix.InMemoryMatrix(data, row_names=row_names, col_names=col_names)
    return M


def get_km_plot_attribute(new_group, name, file_type, geneid):
    """generate new group name and color list for km_plot"""
    if len(name) == 2:
        if file_type == 'gene_expression_file':
            new_name = [single_name.replace(name[0], 'Low Expression')
                        for single_name in new_group]
            new_name = [single_name.replace(name[1], 'High Expression')
                        for single_name in new_name]
            name_update = [single_name.replace(name[0], 'Low Expression')
                           for single_name in name]
            name_update = [single_name.replace(name[1], 'High Expression')
                           for single_name in name_update]
        elif file_type == 'geneset_file':
            new_name = [single_name.replace(
                name[0], 'Low ' + geneids[i] + ' Score')
                for single_name in new_name]
            new_name = [single_name.replace(
                name[1], 'High ' + geneids[i] + ' Score')
                for single_name in new_name]
            name_update = [single_name.replace(
                name[0], 'Low ' + geneids[i] + ' Score')
                for single_name in name]
            name_update = [single_name.replace(
                name[1], 'High ' + geneids[i] + ' Score')
                for single_name in name_update]
        color_command = ('col <- list("' + name_update[0] + '"="#1533AD","'
                         + name_update[1] + '"="#FFB300")')
    elif len(name) > 2:
        if file_type == 'gene_expression_file':
            new_name = [single_name.replace(
                name[0], 'Lowest ' + name[0] + '% Expression')
                for single_name in new_group]
            new_name = [single_name.replace(
                name[-1], 'Highest ' + name[-1] + '% Expression')
                for single_name in new_name]
            name_update = [single_name.replace(
                name[0], 'Lowest ' + name[0] + '% Expression')
                for single_name in name]
            name_update = [single_name.replace(
                name[-1], 'Highest ' + name[-1] + '% Expression')
                for single_name in name_update]
            for k in range(len(name) - 2):
                new_name = [single_name.replace(
                    name[1 + k], 'Middle ' + name[1 + k] + '% Expression')
                    for single_name in new_name]
                name_update = [single_name.replace(
                    name[1 + k], 'Middle ' + name[1 + k] + '% Expression')
                    for single_name in name_update]
        elif file_type == 'geneset_file':
            new_name = [single_name.replace(
                name[0], 'Lowest ' + name[0] + '%' + geneid + ' Score')
                for single_name in new_group]
            new_name = [single_name.replace(
                name[-1], 'Highest ' + name[-1] + '%' + geneid + ' Score')
                for single_name in new_name]
            name_update = [single_name.replace(
                name[0], 'Lowest ' + name[0] + '%' + geneid + ' Score')
                for single_name in name]
            name_update = [single_name.replace(
                name[-1], 'Highest ' + name[-1] + '%' + geneid + ' Score')
                for single_name in name_update]
            for k in range(len(name) - 2):
                new_name = [single_name.replace(
                    name[1 + k], 'Middle ' + name[1 + k] + '% '
                    + geneid + ' Score') for single_name in new_name]
                name_update = [single_name.replace(
                    name[1 + k], 'Middle ' + name[1 + k] + '% '
                    + geneid + ' Score') for single_name in name_update]
        color_list = colorlib.bild_colors(len(name))
        color_command = 'col <- list('
        for k in range(len(color_list)):
            new_color = [hex(int(L * 255))[2:] for L in color_list[k]]
            color_hex = '#' + new_color[0] + new_color[1] + new_color[2]
            color_command = (color_command + '"'
                             + name_update[k] + '"="' + color_hex + '",')
        color_command = color_command[:-1] + ')'
    return new_name, color_command


def main():
    parser = argparse.ArgumentParser(
        description='associate the gene expression with outcomes')
    parser.add_argument(dest='expression_file', type=str, default=None,
                        help='gene expression file or geneset score file')
    parser.add_argument(dest='clinical_data', type=str, default=None)
    parser.add_argument('--outcome', dest='outcome', type=str,
                        help=('format <time_header>,<dead_header>'),
                        default=[], action='append')
    parser.add_argument('--gene', dest='gene', type=str,
                        help=('gene to analyze.  It can appear anywhere '
                              'in the annotations of the expression_file. '
                              'To specify multiple genes, use this parameter '
                              'multiple times.'),
                        default=None, action='append')
    parser.add_argument('--geneset', dest='geneset', type=str,
                        help=('geneset to analyze. '
                              'To specify multiple genes, use this parameter '
                              'multiple times.'),
                        default=[], action='append')
    parser.add_argument('--cutoff', dest='cutoff', type=str,
                        help=('comma-separated list of breakpoints '
                              '(between 0 and 1),e.g. 0.25,0.50'),
                        default=None)
    parser.add_argument('-o', dest='filestem', type=str,
                        help='prefix used to name files.  e.g. "myanalysis".',
                        default=None)
    parser.add_argument('--write_prism', dest='write_prism',
                        action="store_true", default=False,
                        help='write Prism-formatted output')
    parser.add_argument('--plot_km', dest='plot_km', action="store_true",
                        help='write PNG-formatted Kaplan-Meier plot',
                        default=False)
    parser.add_argument('--xlab', dest='xlab', default=False,
                        help='the x label for Kaplan-Meier plot')
    parser.add_argument('--ylab', dest='ylab', default=False,
                        help='the y label for Kaplan-Meier plot')
    args = parser.parse_args()
    input_file = args.expression_file
    assert input_file, ('please specify the path of gene expression data '
                        'file or geneset score file')
    file_type = guess_file_type(input_file)
    clinical_file = args.clinical_data
    assert clinical_file, 'please specify the path of clinical data'
    outcomes = args.outcome
    assert len(outcomes) > 0, 'please specify the time_header and dead_header'
    cutoffs = parse_cutoffs(args.cutoff)
    gene_list = parse_genes(args.gene)
    geneset_list = args.geneset
    if gene_list and geneset_list:
        raise AssertionError('we only accept gene '
                             'or geneset at each time, not both')
    if file_type == 'gene_expression_file' and geneset_list:
        raise AssertionError('gene_expression_file cannot run with geneset')
    if file_type == 'geneset_file' and gene_list:
        raise AssertionError('geneset_file cannot run with gene')
    if file_type == 'gene_expression_file':
        M = arrayio.read(input_file)
    elif file_type == 'geneset_file':
        M = convert_genesetfile2matrix(input_file)
    clinical_data = genesetlib.read_tdf(
        clinical_file, preserve_spaces=True, allow_duplicates=True)
    clinical_dict, name_in_clinical = intersect_clinical_matrix_sample_name(
        M, clinical_data)
    #align the sample order of gene expression file and clinical data
    x = matrixlib.align_cols_to_annot(
        M, name_in_clinical, reorder_MATRIX=True)
    M, colnames = x
    #if given genes or geneset,get of match genes or geneset
    if gene_list:
        x = matrixlib.align_rows_to_annot(M, gene_list, reorder_MATRIX=True)
    if geneset_list:
        x = matrixlib.align_rows_to_annot(M, geneset_list, reorder_MATRIX=True)
    M, rownames = x
    ids = M._row_order
    geneids = M._row_names[ids[0]]
    assert geneids, 'we cannot match any gene or geneset as required'
    headers = ids[:]
    data_all = M.slice()
    #add the gene annotation column to the output_data
    output_data = []
    for i in range(len(geneids)):
        output_data.append([M._row_names[name][i] for name in ids])
    kaplanmeierlib_path = config.kaplanmeierlib
    assert os.path.exists(kaplanmeierlib_path), (
        'can not find the kaplanmeierlib script %s' % kaplanmeierlib_path)
    R = jmath.start_R()
    R('require(splines,quietly=TRUE)')
    R('source("' + config.kaplanmeierlib + '")')
    #generate output headers
    new_cutoffs = [0] + cutoffs + [1]
    name = [''] * (len(new_cutoffs) - 1)
    for j in range(len(new_cutoffs) - 1):
        name[j] = "%d - %d" % (new_cutoffs[j] * 100, new_cutoffs[j + 1] * 100)
    num_samples = ['Num Samples (' + k + ')' for k in name]
    if file_type == 'gene_expression_file':
        ave_expression = ['Average Expression (' + k + ')'for k in name]
    elif file_type == 'geneset_file':
        ave_expression = ['Average gene set score (' + k + ')'for k in name]
    surv90_header = ['90% Survival (' + k + ')' for k in name]
    surv50_header = ['50% Survival (' + k + ')' for k in name]
    newheader = (['p-value'] + num_samples + ave_expression
                 + surv90_header + surv50_header + ['Relationship'])
    #get the time_header, dead_header,time_data,dead_data,
    #and sample_index for each outcome
    all_time_header = []
    all_dead_header = []
    all_sample_index = []
    all_time_data = []
    all_dead_data = []
    for outcome in outcomes:
        x = outcome.split(',')
        assert len(x) == 2, (
            'the outcome format should be <time_header>,<dead_header>')
        time_header, dead_header = x
        time_data = clinical_dict.get(time_header)
        dead_data = clinical_dict.get(dead_header)
        assert time_data, 'there is no column named %s' % time_header
        assert dead_data, 'there is no column named %s' % dead_header
        #only consider the sample who has month and dead information
        sample_index1 = [index for index, item in
                         enumerate(time_data) if len(item) > 0]
        sample_index2 = [index for index, item in
                         enumerate(dead_data) if len(item) > 0]
        sample_index = list(set(sample_index1).intersection(
            set(sample_index2)))
        all_time_header.append(time_header)
        all_dead_header.append(dead_header)
        all_sample_index.append(sample_index)
        all_time_data.append(time_data)
        all_dead_data.append(dead_data)
    #update the output headers
    for k in range(len(outcomes)):
        if len(outcomes) > 1:
            newheader1 = [all_time_header[k] + ' ' + i for i in newheader]
        else:
            newheader1 = newheader
        headers.extend(newheader1)
    #calculate the survival analysis for each gene in each outcome
    all_group_name = [[]] * len(outcomes)
    all_survival = []
    all_dead = []
    all_p_value = [[]] * len(outcomes)
    jmath.R_equals(cutoffs, 'cutoffs')
    for h in range(len(outcomes)):
        survival = [all_time_data[h][i] for i in all_sample_index[h]]
        dead = [all_dead_data[h][i] for i in all_sample_index[h]]
        all_survival.append(survival)
        all_dead.append(dead)
        for i in range(len(geneids)):
            data = data_all[i]
            data_new = [data[j] for j in all_sample_index[h]]
            jmath.R_equals(data_new, 'F')
            R('group <- group.by.value(F, cutoffs)')
            group = R['group']
            avg = [0] * len(name)
            num_group = [0] * len(name)
            group_name = [""] * len(group)
            for k in range(len(name)):
                group_data = [data_new[j] for j in range(len(group))
                              if group[j] == k]
                avg[k] = str(sum(group_data) / len(group_data))
                num_group[k] = str(len(group_data))
            for k in range(len(group)):
                group_name[k] = name[group[k]]
            p_value, surv90, surv50, direction = calc_km(survival,
                                                         dead, group_name)
            if file_type == 'geneset_file':
                direction = direction.replace('expression', 'gene set score')
            single_data = ([str(p_value)] + num_group + avg + surv90
                           + surv50 + [direction])
            output_data[i].extend(single_data)
            all_group_name[h].append(group_name)
            all_p_value[h].append(p_value)
    #generate the prism file and the km plot for each gene in each outcome
    filestem = '' if not args.filestem else args.filestem + '.'
    for h in range(len(outcomes)):
        jmath.R_equals(all_survival[h], 'survival')
        jmath.R_equals(all_dead[h], 'dead')
        for i in range(len(geneids)):
            new_group = ['"' + j + '"' for j in all_group_name[h][i]]
            jmath.R_equals(new_group, 'name')
            if args.write_prism:
                prism_file = str(filestem + all_time_header[h] + '.'
                                 + geneids[i] + '.prism.txt')
                jmath.R_equals('"' + prism_file + '"', 'filename')
                R('write.km.prism.multi(filename,survival, dead, name)')
            if args.plot_km:
                km_plot = (filestem + all_time_header[h] + '.'
                           + geneids[i] + '.km.png')
                new_name, color_command = get_km_plot_attribute(
                    new_group, name, file_type, geneids[i])
                jmath.R_equals(new_name, 'new_name')
                R(color_command)
                jmath.R_equals('"' + km_plot + '"', 'filename')
                jmath.R_equals('"' + geneids[i] + '"', 'title')
                jmath.R_equals('"p_value=' + str(all_p_value[h][i])
                               + '"', 'sub')
                R('xlab<-""')
                R('ylab<-""')
                if args.xlab:
                    jmath.R_equals('"' + args.xlab + '"', 'xlab')
                if args.ylab:
                    jmath.R_equals('"' + args.ylab + '"', 'ylab')
                R('bitmap(file=filename,type="png256")')
                R('plot.km.multi(survival, dead, new_name,'
                  'col=col, main=title, sub=sub,xlab=xlab,ylab=ylab)')
                R('dev.off()')
    #write the output data
    f = sys.stdout
    if args.filestem:
        outfile = args.filestem + '.stats.txt'
        f = file(outfile, 'w')
    print >> f, '\t'.join(headers)
    for i in range(len(geneids)):
        print >> f, '\t'.join(output_data[i])
    f.close()


if __name__ == '__main__':
    main()
