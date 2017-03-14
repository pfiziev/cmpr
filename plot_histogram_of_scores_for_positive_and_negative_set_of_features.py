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
            score = float(buf[4])
            if feature_bins[chrom][start_bin] == SKIP:
                continue
            bins.append((score, feature_bins[chrom][start_bin]))
    return bins


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
plt.figure(figsize=(15, 15))

for in_fname in args.i:
    echo('Processing scores for:', in_fname)
    bins = annotate_bins(in_fname, feature_bins, BIN_SIZE)
    pos = [s for s, l in bins if l == POS]

    neg = [s for s, l in bins if l == NEG]
    to_plot = [pos, neg]
    labels = ['pos ' + str(mean(pos)), 'neg ' + str(mean(neg))]

    skip = [s for s, l in bins if l == SKIP]
    if skip:
        plt.hist(skip, bins=100, normed=1, label='skip' + str(mean(skip)), alpha=0.5)
        to_plot.append(skip)
        labels.append('skip ' + str(mean(skip)))

    plt.hist(to_plot, bins=100, normed=1, label=labels, alpha=0.5)


plt.legend()

plt.savefig(out_fname)









