import sys, math, gzip, re
from utils import *

from itertools import imap

def pearsonr(x, y):
    # Assume len(x) == len(y)
    n = len(x)

    sum_x = float(sum(x))
    sum_y = float(sum(y))

    sum_x_sq = sum(map(lambda x: pow(x, 2), x))
    sum_y_sq = sum(map(lambda x: pow(x, 2), y))

    psum = sum(imap(lambda x, y: x * y, x, y))

    num = psum - (sum_x * sum_y/n)
    den = pow((sum_x_sq - pow(sum_x, 2) / n) * (sum_y_sq - pow(sum_y, 2) / n), 0.5)

    if den == 0:
        return 0

    return num / den


def read_wig(fname, rescale_to_bin_size=None):
    print >>sys.stderr, 'Reading ' + fname

    signal = {}
    chrom = None
    span = None

    with (gzip.open(fname) if fname.endswith('.gz') else open(fname)) as mark_f:
        for line in mark_f:

            if line.startswith('track'):
                continue

            elif line.startswith('fixedStep'):
                chrom = line.split()[1].split('=')[1]
                span = int(re.search(r'span=(\d+)', line).group(1))
                signal[chrom] = []

            else:
                signal[chrom].append(float(line))

    if span is None:
        print 'ERROR: cannot determine span'
        exit(1)

    if rescale_to_bin_size is not None:
        if rescale_to_bin_size % span != 0:
            print 'ERROR: Span must be a multiple of rescale_to_bin_size! ' \
                  'Cannot rescale span:', span, 'to bin_size:', rescale_to_bin_size
            exit(1)

        window = rescale_to_bin_size / span

        print >>sys.stderr, 'Rescaling from:', span, 'to bin_size:', rescale_to_bin_size, 'with window:', window
        signal = dict((chrom, [mean(signal[chrom][i:i + window])
                               for i in xrange(0, len(signal[chrom]), window)])
                      for chrom in signal)

    return signal, span

if len(sys.argv) == 1:
    print 'usage: %s wig-file lag' % __file__
    exit(1)

wig, span = read_wig(sys.argv[1])
lag = int(sys.argv[2])

for l in xrange(lag):
    r = 0
    for chrom in wig:
        if len(wig[chrom]) <= l + 1000:
            continue

        if l > 0:
            signal = wig[chrom][:-l]
        else:
            signal = wig[chrom]
        lagged_signal = wig[chrom][l:]
        r += pearsonr(signal, lagged_signal)

    print 'lag=', l, 'r=', r / len(wig)

