#!/usr/bin/env python

import gzip
from optparse import OptionParser
import os
import sys
import datetime
import math

import re
from utils import mean, std


def echo(*msg):
    print >>sys.stderr, "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), ' '.join(map(str, msg)))


def output_rpkms_for_bed_regions_from_wig(wig_fname, bed_fname):

    wig_data, span = read_wig(wig_fname)

    echo('Calculating RPKM calues')
    log2_rpkms = []
    missing_chroms = set()
    skipped = 0
    scores = {}
    with gzip.open(bed_fname) if bed_fname.endswith('.gz') else open(bed_fname) as in_f:

        for line in in_f:
            buf = line.strip().split('\t')
            state = buf[-1]
            if state not in scores:
                scores[state] = []
            chrom = buf[0]
            if chrom not in wig_data:
                missing_chroms.add(chrom)
                skipped += 1
                continue

            start = int(buf[1]) / span
            end = 1 + (int(buf[2]) - 1) / span
            for bin_idx in xrange(start, end):
                scores[state].append(wig_data[chrom][bin_idx])
    scores = dict((s, mean(scores[s])) for s in scores)
    for s in sorted(scores, key=lambda k: scores[k], reverse=True):
        print '\t'.join([s, str(scores[s])])


def read_wig(fname):
    echo('Reading ' + fname)

    signal = {}
    chrom = None
    span = None

    with (gzip.open(fname) if fname.endswith('.gz') else open(fname)) as mark_f:
        for line in mark_f:

            if line.startswith('track'):
                continue

            elif line.startswith('fixedStep'):
                chrom = re.search(r'chrom=(\w+)', line).group(1)
                _span = int(re.search(r'span=(\w+)', line).group(1))
                if span is not None and span != _span:
                    print 'ERROR: spans don\'t match:', span, _span

                span = _span
                signal[chrom] = []

            else:
                signal[chrom].append(float(line))

    return signal, span



if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-w", "--wig-file", type="string", dest="wig_fname",
                      help="Input wig file with the signal that will be used to calculate the RPKM values "
                           "(can be gzipped bed)", metavar="FILE")

    parser.add_option("-b", "--bed-file", type="string", dest="bed_fname",
                      help="Calculate the FPKM values only for the regions in this bed file.", metavar="FILE")

    (options, args) = parser.parse_args()

    # if no options were given by the user, print help and exit
    if len(sys.argv) == 1:
        parser.print_help()
        exit(0)

    if options.wig_fname:
        output_rpkms_for_bed_regions_from_wig(options.wig_fname,
                                              options.bed_fname)

