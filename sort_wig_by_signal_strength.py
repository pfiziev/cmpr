from utils import *
import sys
import re

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print 'usage: %s wig-file' % __file__
        exit(1)
    data = []
    chrom = None
    step = None
    start = None
    in_file = sys.argv[1]

    with open_file(in_file) as in_f:
        for line in in_f:

            if line.startswith('track'):
                continue

            if line.startswith('fixedStep'):
                chrom = re.search(r'chrom=(\w+)\t', line).group(1)
                start = int(re.search(r'start=(\d+)\t', line).group(1))
                step = int(re.search(r'step=(\d+)\t', line).group(1))
                continue

            data.append((float(line), chrom, start, start + step))
            start += step

    for value, chrom, start, end in sorted(data, reverse=True):
        print '\t'.join(map(str, [chrom, start, end, '.', value, '.']))

    from matplotlib import pyplot as plt
    plt.hist([r[0] for r in data], bins=50)
    plt.savefig(in_file + '.scores_histogram.png')


