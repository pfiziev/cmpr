import argparse
import pybedtools
import sys
from utils import *

genome_sizes = {'hg19': 2897310462,
                'mm9': 2620345972}


def fill_in_missing_states_in_matrix(transition_size, states):
    for s1 in states:
        if s1 not in transition_size:
            transition_size[s1] = {}

        for s2 in states:
            if s2 not in transition_size[s1]:
                transition_size[s1][s2] = 0


def get_transition_sizes(fname):
    states = set()
    transition_size = {}

    echo('Computing transition sizes from:', fname)
    with open_file(fname) as in_f:
        for line in in_f:
            buf = line.strip().split('\t')
            s1, s2 = buf[6].split('-')
            s1 = s1.replace(' <', '')
            s2 = s2.replace('> ', '')
            states.add(s1)
            states.add(s2)
            if s1 not in transition_size:
                transition_size[s1] = {}
            if s2 not in transition_size[s1]:
                transition_size[s1][s2] = 0

            transition_size[s1][s2] += int(buf[2]) - int(buf[1])
    states = sorted(states, key=state_key)

    fill_in_missing_states_in_matrix(transition_size, states)

    return transition_size, states

def get_overlaps(fgr_bed, chrom_compare_bed, states):
    overlaps = {}
    for f in chrom_compare_bed.intersect(b=fgr_bed).features():
        s1, s2 = f[6].split('-')
        s1 = s1.replace(' <', '')
        s2 = s2.replace('> ', '')

        if s1 not in overlaps:
            overlaps[s1] = {}
        if s2 not in overlaps[s1]:
            overlaps[s1][s2] = 0

        overlaps[s1][s2] += int(f[2]) - int(f[1])

    fill_in_missing_states_in_matrix(overlaps, states)
    return overlaps

def get_total_coverage(fname):
    coverage = 0
    echo('Computing coverage in:', fname)
    with open_file(fname) as in_f:
        for line in in_f:
            buf = line.strip().split('\t')
            coverage += int(buf[2]) - int(buf[1])
    echo('Coverage:', coverage)
    return coverage

PSEUDO_COUNT = 1

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-a', help='chrom_compare output (bed file)')
    parser.add_argument('-b', help='bed file with regions to overlap')
    parser.add_argument('-B', help='Background for -b (optional)')

    parser.add_argument('-g', '--genome', choices=['hg19', 'mm9'])

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)

    chrom_compare_outfname = args.a

    transition_size, states = get_transition_sizes(chrom_compare_outfname)
    chrom_compare_bed = pybedtools.BedTool(chrom_compare_outfname)

    fgr_bed = pybedtools.BedTool(args.b)
    fgr_size = get_total_coverage(args.b)
    bgr_bed = None
    bgr_size = None

    if args.B:
        echo('Using background from:', args.B)
        bgr_bed = pybedtools.BedTool(args.B)
        bgr_size = get_total_coverage(args.B)
    else:
        bgr_size = genome_sizes[args.genome]

    echo('Background size:', bgr_size)

    fgr_overlaps = get_overlaps(fgr_bed, chrom_compare_bed, states)
    if bgr_bed is not None:
        bgr_overlaps = get_overlaps(bgr_bed, chrom_compare_bed, states)

    print '\t'.join(['State'] + states)
    for s1 in states:
        enrichments = []
        for s2 in states:
            if bgr_bed is None:
                enr = (fgr_overlaps[s1][s2] + PSEUDO_COUNT) / float(PSEUDO_COUNT + transition_size[s1][s2] * fgr_size / bgr_size)
            else:
                enr = (fgr_overlaps[s1][s2] + PSEUDO_COUNT)/ float(PSEUDO_COUNT + bgr_overlaps[s1][s2] * fgr_size / bgr_size)
            enrichments.append(enr)
        print '\t'.join([s1] + ['%.2lf' % e for e in enrichments])



