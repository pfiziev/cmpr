import argparse
from itertools import izip
import os
import pprint
import sys
import math
import pybedtools
from utils import echo, mean, open_file, standardize
import matplotlib
matplotlib.use('Agg')

__author__ = 'Fiziev'


def read_info(fnames):
    info = {}
    echo('Reading info from:', fnames)
    for fname in fnames:
        group = os.path.split(fname)[1].split('.')[0]
        for l in open_file(fname):
            dataset = l.strip().split()[0]
            info[dataset] = group

    return info


def read_gene_expression(gene_info_fname, expression_fname, group_info, TSS_WINDOW):

    gene_info = []
    expression = {}

    with open_file(expression_fname) as in_f:

        celltypes = in_f.readline().strip().split()[1:]

        for line in in_f:
            buf = line.strip().split()
            gene_id = buf[0]
            expression[gene_id] = {}
            cur_gene_expression = map(lambda v: math.log(float(v) + 1, 2), buf[1:])

            for ct, expr in izip(celltypes, cur_gene_expression):

                if ct not in group_info:
                    continue

                group_id = group_info[ct]

                if group_id not in expression[gene_id]:
                    expression[gene_id][group_id] = []

                expression[gene_id][group_id].append(expr)

            for group_id in expression[gene_id]:
                expression[gene_id][group_id] = mean(expression[gene_id][group_id])

    with open_file(gene_info_fname) as in_f:

        for l in in_f:
            buf = l.strip().split()
            gene_id = buf[0]

            if gene_id not in expression:
                continue

            if buf[1].startswith('GL') or buf[1] == 'M':
                continue

            chrom = 'chr' + buf[1]

            tss = int(buf[3] if buf[4] == '-1' else buf[2])

            gene_info.append('\t'.join(map(str, [chrom, max(0, tss - TSS_WINDOW), tss + TSS_WINDOW + 1, gene_id])))

    gene_info = pybedtools.BedTool('\n'.join(gene_info), from_string=True)
    # print_diff_genes_stats(expression)
    return gene_info, expression


def print_diff_genes_stats(expression):
    groups = sorted(expression.values()[0])
    print 'all genes:', len(expression)
    for i, g1 in enumerate(groups):
        for g2 in groups[i + 1:]:
            print '\t'.join(map(str, [g1, g2,
                                      sum(is_differential(gene_id, g1, g2, expression, diff_level=1) for gene_id in expression),
                                      sum(is_differential(gene_id, g1, g2, expression, diff_level=2) for gene_id in expression),
                                      sum(is_differential(gene_id, g1, g2, expression, diff_level=3) for gene_id in expression),
                                      sum(is_differential(gene_id, g1, g2, expression, diff_level=4) for gene_id in expression),
                                      sum(is_differential(gene_id, g1, g2, expression, diff_level=5) for gene_id in expression),
                                      sum(is_differential(gene_id, g1, g2, expression, diff_level=6) for gene_id in expression),
                                      sum(is_differential(gene_id, g1, g2, expression, diff_level=7) for gene_id in expression)
                                      ]))



def is_differential(gene_id, group_A, group_B, gene_expression, diff_level):
    return abs(gene_expression[gene_id][group_A] - gene_expression[gene_id][group_B]) >= diff_level


FDR = 'fdr'
SCORE = 'score'

def annotate_bins(fname, gene_info, gene_expression, diff_level):

    bins = []

    group_A, group_B = os.path.split(fname)[1].split('.')[0].split('_vs_')
    groups = set(gene_expression.values()[0])

    if group_A not in groups or group_B not in groups:
        echo('No expression info for:', fname)
        return []

    for f in pybedtools.BedTool(fname).intersect(b=gene_info, wo=True):

        score = float(f[4])
        fdr = float(f[8])
        gene_id = f[12]

        bins.append((score, fdr, is_differential(gene_id, group_A, group_B, gene_expression, diff_level)))

    return bins



parser = argparse.ArgumentParser()

parser.add_argument('-i', nargs='+', help='roadmap info files')
parser.add_argument('-q', nargs='+', help='q-value scores')
parser.add_argument('-e', help='expression matrix')
parser.add_argument('-g', help='gene info')
parser.add_argument('-w', help='TSS window size', default=2000, type=int)
parser.add_argument('--diff-level', help='delta log2 rpkm', default=1, type=float)

parser.add_argument('-o', '--output', help='output prefix')

args = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    exit(1)

echo('Window:', args.w)
group_info = read_info(args.i)
gene_info, gene_expression = read_gene_expression(args.g, args.e, group_info, args.w)
out_fname = args.output
diff_level = args.diff_level

RAW_SCORES = 'raw_scores'
RAW_SCORES_SMOOTH = 'raw_scores_smooth'
q_values_info = {RAW_SCORES: {},
                 RAW_SCORES_SMOOTH: {}}
for fname in args.q:
    group_id, method = os.path.split(fname)[1].split('.')[:2]
    if method not in q_values_info:
        q_values_info[method] = {}
    q_values_info[method][group_id] = fname
    if 'smooth' in fname:
        q_values_info[RAW_SCORES_SMOOTH][group_id] = fname
    else:
        q_values_info[RAW_SCORES][group_id] = fname

from matplotlib import pyplot as plt
plt.figure(0, figsize=(10, 10))
plt.figure(1, figsize=(10, 10))

method_to_label = {'shuffle_segmentations_use_original_combo_for_background_smooth': 'seg+real/sm',
                              'shuffle_segmentations_use_original_combo_for_background': 'seg+real',
                              RAW_SCORES: 'raw',
                              RAW_SCORES_SMOOTH: 'raw/sm',
                              'shuffle_marks': 'marks',
                              'shuffle_marks_smooth': 'marks/sm',
                              'shuffle_segmentations': 'seg',
                              'shuffle_segmentations_smooth': 'seg/sm',
                              }

for method in sorted(q_values_info):
    echo('Processing scores for:', method)
    bins = []
    for group_id in q_values_info[method]:
        fname = q_values_info[method][group_id]
        echo(fname)
        bins.extend(annotate_bins(fname, gene_info, gene_expression, diff_level))

    echo('Computing ROC curve')
    TPR = []
    FPR = []
    precision = []

    total_positive = sum(is_diff for _, _, is_diff in bins)
    total_negative = len(bins) - total_positive

    echo('Total positive:', total_positive)
    echo('Total negative:', total_negative)

    fp_so_far = 0
    tp_so_far = 0
    auc = 0
    auc_PR = 0
    for _, _, is_diff in sorted(bins,
                                key=lambda (score, fdr, _): score if method in [RAW_SCORES, RAW_SCORES_SMOOTH] else (fdr, -score),
                                reverse=(method in [RAW_SCORES, RAW_SCORES_SMOOTH])):
        if is_diff:
            tp_so_far += 1
        else:
            fp_so_far += 1

        TPR.append(tp_so_far / float(total_positive))
        FPR.append(fp_so_far / float(total_negative))
        precision.append(tp_so_far / float(tp_so_far + fp_so_far))

        if len(FPR) > 1:
            auc += (FPR[-1] - FPR[-2]) * (TPR[-2] + TPR[-1]) / 2.
            auc_PR += (TPR[-1] - TPR[-2]) * (precision[-2] + precision[-1]) / 2.

    plt.figure(0)
    plt.plot(FPR, TPR, label=method_to_label[method] + ', %.4lf' % auc)

    plt.figure(1)
    plt.plot(TPR, precision, label=method_to_label[method] + ', %.4lf' % auc_PR)


plt.figure(0)
plt.plot([0, 1], [0, 1], color='black')

plt.xlabel('FPR')
plt.ylabel('TPR')

plt.legend(loc='lower right')

plt.savefig(out_fname + '.ROC.png')

plt.figure(1)

plt.xlabel('Recall')
plt.ylabel('Precision')

plt.legend(loc='upper right')

plt.savefig(out_fname + '.PR.png')










