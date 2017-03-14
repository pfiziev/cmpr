import argparse
from itertools import izip
import os
import pprint
import random
import sys
import math
import pybedtools
from utils import echo, mean, open_file, standardize
import matplotlib
matplotlib.use('Agg')

__author__ = 'Fiziev'

FDR = 'fdr'
SCORE = 'score'



POS = 1
NEG = -1
SKIP = 0

def read_features(pos_fname, neg_fname, chrom_lengths, bin_size):
    echo('Reading features:', pos_fname, neg_fname)
    features = dict((chrom, [SKIP] * chrom_lengths[chrom]) for chrom in chrom_lengths)

    with open_file(pos_fname) as in_f:
        for l in in_f:
            chrom, start, end = l.strip().split('\t')[:3]
            start_bin = int(start) / bin_size
            end_bin = int(end) / bin_size

            if chrom not in features:
                features[chrom] = []

            for bin_i in xrange(start_bin, end_bin + 1):
                if bin_i >= len(features[chrom]):
                    continue
                features[chrom][bin_i] = POS

    if neg_fname is not None:
        with open_file(neg_fname) as in_f:
            for l in in_f:
                chrom, start, end = l.strip().split('\t')[:3]
                start_bin = int(start) / bin_size
                end_bin = int(end) / bin_size

                if chrom not in features:
                    features[chrom] = []

                for bin_i in xrange(start_bin, end_bin + 1):
                    if bin_i >= len(features[chrom]):
                        continue

                    features[chrom][bin_i] = NEG
    else:
        for chrom in features:
            for bin_i in xrange(len(features[chrom])):
                if features[chrom][bin_i] != POS:
                    features[chrom][bin_i] = NEG
    return features


def get_chrom_lengths(fname, bin_size):
    echo('Computing chromosome lengths:', fname)
    chrom_lengths = {}
    with open_file(fname) as in_f:
        for l in in_f:
            chrom, start, end = l.strip().split()[:3]
            if chrom not in chrom_lengths:
                chrom_lengths[chrom] = 0

            chrom_lengths[chrom] = max(chrom_lengths[chrom], int(end) / bin_size)
    return chrom_lengths


def annotate_bins(in_fname, feature_bins, bin_size):
    bins = []
    with open_file(in_fname) as in_f:
        for l in in_f:
            buf = l.strip().split('\t')
            chrom = buf[0]
            start_bin = int(buf[1]) / bin_size
            # end_bin = int(buf[2]) / bin_size
            # score = abs(float(buf[4]))
            score = float(buf[4])

            # if score == 0.0:
            #     continue

            if feature_bins[chrom][start_bin] == SKIP:
                continue
            bins.append((chrom, start_bin, score, feature_bins[chrom][start_bin]))
    return [(score, f) for _, _, score, f in sorted(bins)]


parser = argparse.ArgumentParser()

parser.add_argument('-i', nargs='+', help='score files')
parser.add_argument('-p', help='positive features')
parser.add_argument('-n', help='negative features')

parser.add_argument('-b', type=int, default=200, help='bin size')

parser.add_argument('-o', '--output', help='output prefix')

args = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    exit(1)

BIN_SIZE = args.b
out_fname = args.output

chrom_lengths = get_chrom_lengths(args.i[0], BIN_SIZE)
print chrom_lengths
feature_bins = read_features(args.p, args.n, chrom_lengths, BIN_SIZE)

from matplotlib import pyplot as plt
plt.figure(0, figsize=(25, 25))
plt.figure(1, figsize=(25, 25))

for in_fname in args.i:
    echo('Processing scores for:', in_fname)
    bins = annotate_bins(in_fname, feature_bins, BIN_SIZE)
    random.seed(10000000000)
    bins = random.sample(bins, min(3000000, len(bins)))
    # random.shuffle(bins)


    echo('Computing ROC curve')
    TPR = []
    FPR = []
    precision = []

    total_positive = sum(f == POS for _, f in bins)
    total_negative = len(bins) - total_positive

    echo('Total positive:', total_positive)
    echo('Total negative:', total_negative)

    fp_so_far = 0
    tp_so_far = 0
    auc = 0
    auc_PR = 0
    for _, f in sorted(bins, reverse=True, key=lambda (s, l): s):
        if f == POS:
            tp_so_far += 1
        else:
            fp_so_far += 1
        tpr = tp_so_far / float(total_positive)
        fpr = fp_so_far / float(total_negative)
        prec = tp_so_far / float(tp_so_far + fp_so_far)
        TPR.append(tp_so_far / float(total_positive))
        FPR.append(fp_so_far / float(total_negative))
        precision.append(tp_so_far / float(tp_so_far + fp_so_far))

        if len(FPR) > 1:
            auc += (FPR[-1] - FPR[-2]) * (TPR[-2] + TPR[-1]) / 2.
            auc_PR += (TPR[-1] - TPR[-2]) * (precision[-2] + precision[-1]) / 2.

    import numpy as np
    x = np.linspace(0, 1, num=1000)

    # plt.figure(0)
    # plt.plot(x, np.interp(x, FPR, TPR), label=in_fname + ', %.4lf' % auc)
    #
    # plt.figure(1)
    # plt.plot(x, np.interp(x, TPR, precision), label=in_fname + ', %.4lf' % auc_PR)

    plt.figure(0)
    plt.plot(FPR, TPR, label=in_fname + ', %.4lf' % auc)

    plt.figure(1)
    plt.plot(TPR, precision, label=in_fname + ', %.4lf' % auc_PR)


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










