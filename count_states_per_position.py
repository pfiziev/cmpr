import argparse
from utils import *
import operator
from metric import *

BIN_SIZE = 200

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-a', nargs='+', help='segmentation files of group A')

    parser.add_argument('-g', help='group file')

    parser.add_argument('-o', help='output prefix')

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)

    if args.g:
        seg_fnames = [os.path.join(os.path.split(args.g)[0], l.strip()) for l in open(args.g)]
    else:
        seg_fnames = args.a

    segmentations, states = read_segmentations(seg_fnames)
    out_prefix = args.o

    chromosomes = sorted(reduce(operator.and_, [set(s.keys())
                                                for s in segmentations.values()]),
                         key=lambda c: int(c.replace('chr', '')) if re.search(r'^chr(\d+)$', c) else 100)

    # print chromosomes
    segmentations = filter_chroms(segmentations, chromosomes)

    wig_files = dict((s, gzip.open(out_prefix + '.' + s + '.state_counts.wig.gz', 'w')) for s in states)

    for s in wig_files:
        title = os.path.split(out_prefix)[1] + '.' + s
        wig_files[s].write('track type=wiggle_0 name="%s" description="%s"\n' % (title, title))

    for chrom in chromosomes:
        echo('Processing:', chrom)

        for s in wig_files:
            wig_files[s].write('fixedStep\tchrom=%s\tstart=0\tstep=%d\tspan=%d\n' % (chrom, BIN_SIZE, BIN_SIZE))

        for i, bins in enumerate(izip(*[iterstate(seg[chrom], BIN_SIZE) for seg in segmentations.values()])):
            state_counts = dict((s, 0) for s in states)
            for s in bins:
                state_counts[s] += 1
            for s in state_counts:
                wig_files[s].write('%d\n' % state_counts[s])

    for s in wig_files:
        wig_files[s].close()

    echo('Done')

