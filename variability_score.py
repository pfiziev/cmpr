import argparse
from itertools import izip
import os
import sys
import math
import operator
from utils import *
from metric import *


def get_variability_score(chrom, segmentations, metric, BIN_SIZE):
    score = []
    datasets = sorted(segmentations)

    # print chrom
    for bins in izip(*[iterstate(segmentations[d][chrom], BIN_SIZE) for d in datasets]):

        mean_probs = get_prob_vector(bins, datasets, metric)
        # variability = sum((symmetric_KL_divergence(mean_probs)))
        variability = sum((KL(P=metric[dataset][state],
                              Q=mean_probs) ** 2
                           for ((_, state), dataset) in izip(bins, datasets))) / len(datasets)

        score.append(variability)
        # print 'A:', a_bins
        # print 'B:', b_bins
        # print

    return smooth(score)

def store_wig(chrom, wig_data, out_f, span):
    out_f.write('fixedStep\tchrom=%s\tstart=0\tstep=%d\tspan=%d\n' % (chrom, span, span))
    for v in wig_data:
        out_f.write('%.2lf\n' % v)


parser = argparse.ArgumentParser()

parser.add_argument('-a', nargs='+', help='segmentation files')
parser.add_argument('--bin-size', type=int, help='Bin size, default: %(default)s', default=200)

parser.add_argument('-o', help='output file')


args = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    exit(1)

BIN_SIZE = args.bin_size


segmentations, states = read_segmentations(args.a)

chromosomes = sorted(reduce(operator.and_, [set(s.keys()) for s in segmentations.values()]))
# print chromosomes
segmentations = filter_chroms(segmentations, chromosomes)

metric = learn_metric(segmentations, states, BIN_SIZE)
#print_metric(metric_A)

out_fname = args.o

with open_file(out_fname + '.wig.gz', 'w') as out_f:
    title = os.path.split(out_fname)[1]
    out_f.write('track type=wiggle_0 name="%s" description="%s"\n' % (title, title))
    for chrom in chromosomes:
        if not all(chrom in segmentations[d] for d in segmentations):
            echo(chrom, 'is not found in all segmentations. Skipping..')
            continue
        store_wig(chrom,
                  get_variability_score(chrom, segmentations, metric, BIN_SIZE),
                  out_f=out_f,
                  span=BIN_SIZE)


