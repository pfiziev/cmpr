import random
import sys
from utils import open_file

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print 'usage: %s segmentation.bed [chrom]' % __file__
        exit(1)

    filter_chrom = None
    if len(sys.argv) == 3:
        filter_chrom = sys.argv[2]

    state_lengths = {}
    bin_size = 200
    segmentation = []
    with open_file(sys.argv[1]) as in_f:
        for line in in_f:
            chrom, start, end, state = line.strip().split('\t')

            if filter_chrom:
                if chrom != filter_chrom:
                    continue

            if state not in state_lengths:
                state_lengths[state] = 0
            reg_length = int(end) - int(start)
            segmentation.extend([state] * (reg_length / bin_size))

    random.shuffle(segmentation)

    for bin_no, state in enumerate(segmentation):
        print '\t'.join(map(str, ['chrZ' if filter_chrom is None else filter_chrom,
                                  bin_no * bin_size,
                                  (bin_no + 1) * bin_size,
                                  state]))

