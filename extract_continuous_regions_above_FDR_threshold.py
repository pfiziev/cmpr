import argparse
import sys
from utils import open_file

__author__ = 'Fiziev'


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', help='scores')
    parser.add_argument('-f', '--fdr-threshold', type=float, help='FDR threshold: %(default)s', default=0.01)

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)

    data = {}
    with open_file(args.i) as in_f:
        for line in in_f:
            buf = line.strip().split()
            chrom, start, end, score, state_transition, fdr = buf[0], buf[1], buf[2], buf[4], buf[6], buf[8]
            start = int(start)
            end = int(end)
            score = float(score)
            fdr = float(fdr)

            if fdr <= args.fdr_threshold:
                if chrom not in data:
                    data[chrom] = []

                data[chrom].append((start, end, state_transition, score))

    for chrom in sorted(data):
        regs = sorted(data[chrom])
        locus_start = None
        locus_score = 0
        locus_length = 0

        for i in xrange(len(regs)):

            start, end, state_transition, score = regs[i]

            if locus_start is None:
                locus_start = start
                locus_score = 0
                locus_length = 0

            locus_score += score
            locus_length += 1

            if i == len(regs) - 1 or regs[i + 1][0] != regs[i][1] or regs[i + 1][2] != state_transition:
                print '\t'.join(map(str,
                                    [chrom,
                                     locus_start,
                                     end,
                                     '.',
                                     locus_score / locus_length,
                                     '.',
                                     state_transition]))
                locus_start = None





