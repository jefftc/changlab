#!/usr/bin/env python
from genomicode import filelib, jmath, genesetlib, config, matrixlib
import arrayio
import os
import argparse
import numpy
import sys


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


def main():
    parser = argparse.ArgumentParser(
        description='associate the gene expression with outcomes')
    parser.add_argument(dest='expression_file', type=str, default=None)
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
    args = parser.parse_args()
    input_file = args.expression_file
    assert input_file, 'please specify the path of gene expression data file'
    clinical_file = args.clinical_data
    assert clinical_file, 'please specify the path of clinical data'
    outcomes = args.outcome
    assert len(outcomes) > 0, 'please specify the time_header and dead_header'
    cutoffs = args.cutoff
    if cutoffs:
        cutoffs = cutoffs.split(',')
        cutoffs = [float(i) for i in cutoffs]
        cutoffs.sort()
    if not cutoffs:
        cutoffs = [0.5]
    for cutoff in cutoffs:
        assert cutoff > 0 and cutoff < 1.0, 'cutoff should be between 0 and 1'
    genes = args.gene
    if genes:
        gene_list = []
        for i in genes:
            gene_list.extend(i.split(','))
    M = arrayio.read(input_file)
    clinical_data = genesetlib.read_tdf(
        clinical_file, preserve_spaces=True, allow_duplicates=True)
    sample_name = M._col_names['_SAMPLE_NAME']
    name_in_clinical = None
    clinical_dict = dict()
    for g in clinical_data:
        clinical_dict[g[0]] = g[2]
    for g in clinical_dict:
        n = len(list(set(sample_name).intersection(set(clinical_dict[g]))))
        if n > 0:
            name_in_clinical = clinical_dict[g]
    assert name_in_clinical, ('cannot match the sample name in '
                              'clinical data and gene expression data')
    #align the sample order of gene expression file and clinical data
    x = matrixlib.align_cols_to_annot(
        M, name_in_clinical, reorder_MATRIX=True)
    M, colnames = x
    #if given genes,get the row index of match genes
    if genes:
        row_index = []
        for gene in gene_list:
            for column_name in M._row_order:
                if gene in M.row_names(column_name):
                    index_list = [index for index, item in
                                  enumerate(M.row_names(column_name))
                                  if item == gene]
                    row_index.extend(index_list)
        assert row_index, ('we cannot find the genes you given '
                           'in the expression data')
        M = M.matrix(row_index, None)
    ids = M._row_order
    geneid = M._row_names[ids[0]]
    headers = ids[:]
    data_all = M.slice()
    output_data = []
    for i in range(len(geneid)):
        output_data.append([M._row_names[name][i] for name in ids])
    kaplanmeierlib_path = config.kaplanmeierlib
    assert os.path.exists(kaplanmeierlib_path), (
        'can not find the kaplanmeierlib script %s' % kaplanmeierlib_path)
    R = jmath.start_R()
    R('require(splines,quietly=TRUE)')
    R('source("' + config.kaplanmeierlib + '")')
    for outcome in outcomes:
        time_data = None
        dead_data = None
        x = outcome.split(',')
        assert len(x) == 2, (
            'the outcome format should be <time_header>,<dead_header>')
        time_header, dead_header = x
        for g in clinical_dict:
            if g == time_header:
                time_data = clinical_dict[g]
            if g == dead_header:
                dead_data = clinical_dict[g]
        assert time_data, 'there is no column named %s' % time_header
        assert dead_data, 'there is no column named %s' % dead_header
        #only consider the sample who has month and dead information
        sample_index1 = [index for index, item in
                         enumerate(time_data) if len(item) > 0]
        sample_index2 = [index for index, item in
                         enumerate(dead_data) if len(item) > 0]
        sample_index = list(set(sample_index1).intersection(
            set(sample_index2)))
        new_cutoffs = [0]
        new_cutoffs.extend(cutoffs)
        new_cutoffs.append(1)
        name = [''] * (len(new_cutoffs) - 1)
        for j in range(len(new_cutoffs) - 1):
            name[j] = str(str(int(new_cutoffs[j] * 100)) + '-' +
                          str(int(new_cutoffs[j + 1] * 100)) + '%')
        newheader = ['p-value']
        num_samples = ['Num Samples (' + k + ')' for k in name]
        ave_expression = ['Average Expression (' + k + ')'for k in name]
        surv90_header = ['90% Survival (' + k + ')' for k in name]
        surv50_header = ['50% Survival (' + k + ')' for k in name]
        newheader.extend(num_samples)
        newheader.extend(ave_expression)
        newheader.extend(surv90_header)
        newheader.extend(surv50_header)
        newheader.append('Relationship')
        if len(outcomes) > 1:
            newheader = [time_header + ' ' + i for i in newheader]
        headers.extend(newheader)
        survival = [time_data[i] for i in sample_index]
        dead = [dead_data[i] for i in sample_index]
        jmath.R_equals(cutoffs, 'cutoffs')
        for i in range(len(geneid)):
            data = data_all[i]
            data_new = [data[j] for j in sample_index]
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
            single_data = [str(p_value)]
            single_data.extend(num_group)
            single_data.extend(avg)
            single_data.extend(surv90)
            single_data.extend(surv50)
            single_data.append(direction)
            output_data[i].extend(single_data)
            filestem = '' if not args.filestem else args.filestem + '.'
            new_group = ['"' + j + '"' for j in group_name]
            jmath.R_equals(new_group, 'name')
            if args.write_prism:
                prism_file = str(filestem + time_header + '.' + geneid[i] +
                                 '.prism.txt')
                jmath.R_equals('"' + prism_file + '"', 'filename')
                R('write.km.prism.multi(filename,survival, dead, name)')
            if args.plot_km:
                km_plot = filestem + time_header + '.' + geneid[i] + '.km.png'
                jmath.R_equals('"' + km_plot + '"', 'filename')
                R('col <- list("' + name[0] + '"="#FF0000")')
                R('bitmap(file=filename,type="png256")')
                R('plot.km.multi(survival, dead, name, col=col)')
                R('dev.off()')
    f = sys.stdout
    if args.filestem:
        outfile = args.filestem + '.stats.txt'
        f = file(outfile, 'w')
    print >> f, '\t'.join(headers)
    for i in range(len(geneid)):
        print >> f, '\t'.join(output_data[i])
    f.close()


if __name__ == '__main__':
    main()
