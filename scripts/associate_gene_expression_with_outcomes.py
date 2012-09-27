#!/usr/bin/env python
from genomicode import filelib, jmath, genesetlib, config, matrixlib
import arrayio
import os
import argparse
import numpy
import sys


def main():
    parser = argparse.ArgumentParser(
        description='associate the gene expression with outcomes')
    parser.add_argument(dest='expression_file', type=str,
                        help='specify the gene expression data file',
                        default=None)
    parser.add_argument(dest='clinical_data', type=str,
                        help='specify the clinical data file', default=None)
    parser.add_argument('--outcome', dest='outcome', type=str,
                        help=('specify time_header and dead_header '
                              'in the format<time_header>,<dead_header>'),
                        default=[], action='append')
    parser.add_argument('--genes', dest='genes', type=str,
                        help='specify genes for analysis', default=None,
                        action='append')
    parser.add_argument('-o', dest='outpath', type=str,
                        help='specify the outpath', default=None)
    parser.add_argument('--prism_file', dest='prism_file', type=str,
                        help='specify the prism file name', default=None)
    parser.add_argument('--cutoff', dest='cutoff', type=float,
                        help=('specify the cutoff(between 0 and 1) '
                              'to calculate'),
                        default=[], action='append')
    parser.add_argument('--km_plot', dest='km_plot', type=str,
                        help='specify the file name for km plot', default=None)
    args = parser.parse_args()
    input_file = args.expression_file
    assert input_file, 'please specify the path of gene expression data file'
    clinical_file = args.clinical_data
    assert clinical_file, 'please specify the path of clinical data'
    outcomes = args.outcome
    assert len(outcomes) > 0, 'please specify the time_header and dead_header'
    cutoffs = args.cutoff
    if not cutoffs:
        cutoffs = [0.5]
    for cutoff in cutoffs:
        assert cutoff > 0 and cutoff < 1.0, 'cutoff should be between 0 and 1'
    genes = args.genes
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
    R('ow<-options("warn")')
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
        jmath.R_equals(survival, 'survival')
        dead = [dead_data[i] for i in sample_index]
        jmath.R_equals(dead, 'dead')
        jmath.R_equals(cutoffs, 'cutoffs')
        for i in range(len(geneid)):
            data = data_all[i]
            data_new = [data[j] for j in sample_index]
            jmath.R_equals(data_new, 'F')
            R('group <- group.by.value(F, cutoffs)')
            group = R['group']
            R('name <- rep("", length(survival))')
            avg = [0] * len(name)
            num_group = [0] * len(name)
            for k in range(len(name)):
                jmath.R_equals('"' + name[k] + '"',
                               'name[group ==' + str(k) + ']')
                group_data = [data_new[j] for j in range(len(group))
                              if group[j] == k]
                avg[k] = str(sum(group_data) / len(group_data))
                num_group[k] = str(len(group_data))
            R('x <- calc.km.multi(survival, dead, name)')
            R('options(ow)')
            c = R['x']
            p_value = c.rx2('p.value')[0]
            surv90 = [''] * len(name)
            surv50 = [''] * len(name)
            for k in range(len(name)):
                surv90[k] = c.rx2('surv').rx2(name[k]).rx2('surv.90')[0]
                surv50[k] = c.rx2('surv').rx2(name[k]).rx2('surv.50')[0]
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
                direction = ('Low expression has worse survival.')
            elif med_high < med_low:
                direction = ('High expression has worse survival.')
            if p_value >= 0.05:
                direction = ''                                   
            surv90 = ['' if (str(k) == 'NA') else str(k) for k in surv90]
            surv50 = ['' if (str(k) == 'NA') else str(k) for k in surv50]
            single_data = [str(p_value)]
            single_data.extend(num_group)
            single_data.extend(avg)
            single_data.extend(surv90)
            single_data.extend(surv50)
            single_data.append(direction)
            output_data[i].extend(single_data)
            if args.prism_file:
                assert len(geneid) == 1, ('multiple genes match and cannot '
                                          'write the prism file')
                jmath.R_equals('"' + args.prism_file + '"', 'filename')
                R('write.km.prism.multi(filename,survival, dead, name)')
            if args.km_plot:
                assert len(geneid) == 1, ('multiple genes match and cannot '
                                          'plot the km file')
                jmath.R_equals('"' + args.km_plot + '"', 'filename')
                R('col <- list("' + name[0] + '"="#FF0000")')
                R('pdf(filename)')
                R('plot.km.multi(survival, dead, name, col=col)')
                R('dev.off()')
    f = sys.stdout
    if args.outpath:
        f = file(args.outpath, 'w')
    print >> f, '\t'.join(headers)
    for i in range(len(geneid)):
        print >> f, '\t'.join(output_data[i])
    f.close()


if __name__ == '__main__':
    main()
