import argparse
import sys
import operator
from utils import *
from metric import *

def compute_distances(segmentations, states, BIN_SIZE):
    distances = dict((d1, dict((d2, None) for d2 in segmentations)) for d1 in segmentations)
    for d1 in segmentations:
        for d2 in segmentations:
            if distances[d1][d2] is None:
                distances[d1][d2] = distances[d2][d1] = log_enrichment_distance(segmentations[d1],
                                                                                segmentations[d2],
                                                                                states,
                                                                                BIN_SIZE)

                echo('Computing:', d1, 'vs', d2, 'Distance:', distances[d1][d2])

    return distances

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-a', nargs='+', help='segmentation files of group A')
    parser.add_argument('--bin-size', type=int, help='Bin size, default: %(default)s', default=200)

    parser.add_argument('-o', help='output prefix')
    parser.add_argument('-l', help='labels file')

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)

    BIN_SIZE = args.bin_size

    segmentations, states = read_segmentations(args.a)

    chromosomes = sorted(reduce(operator.and_, [set(s.keys()) for s in segmentations.values()]))
    # print chromosomes
    segmentations = filter_chroms(segmentations, chromosomes)

    distances = compute_distances(segmentations, states, BIN_SIZE)
    if args.l:
        labels = dict(l.strip().split() for l in open(args.l))
    else:
        labels = None

    def label(fname):
        if labels is not None:
            return labels[os.path.split(fname)[1]]
        else:
            return re.sub(r'_(\d+)_segments.bed(.gz)*', '', os.path.split(fname)[1])

    with open_file(args.o, 'w') as out_f:
        out_f.write('\t'.join(['Dataset'] + [label(fn) for fn in args.a]) + '\n')
        for d1 in args.a:
            out_f.write('\t'.join([label(d1)] + map(str, [distances[d1][d2] for d2 in args.a])) + '\n')


    # print_metric(metric_A)
