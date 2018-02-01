import sys, argparse
from utils import *


def count_freq(data, states):
    h = dict((s1, dict((s2, 0) for s2 in states)) for s1 in states)
    for _, s1, s2 in data:
        h[s1][s2] += 1
    return h

parser = argparse.ArgumentParser()

parser.add_argument('-i', nargs='+', help='score file (chrom_compare output: .bed.gz)')
parser.add_argument('-t', type=float, help='threshold (on FDR by default)')
parser.add_argument('-p', action='store_true', help='use top %% instead of FDR cutoff')
parser.add_argument('-s', action='store_true', help='use scores instead of FDR for cutoff')

args = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    exit(1)


def compute_state_change_fractions(in_fname, symmetric=False):
    data = []
    states = set()
    with open_file(in_fname) as in_f:
        for l in in_f:
            buf = l.strip().split()
            score = abs(float(buf[4]))
            s1, s2 = buf[6].split('-')
            qvalue = float(buf[8])
            states.add(s1)
            states.add(s2)

            data.append((score if args.s else qvalue, s1, s2))
            if symmetric:
                data.append((score if args.s else qvalue, s2, s1))
    states = sorted(states, key=state_key)
    global_stats = count_freq(data, states)
    n_total = float(len(data))
    if args.p:
        n = len(data)
        data = sorted(data, reverse=True)[:int(n * args.t)]
    else:
        data = filter(lambda (s, x, y): (s >= args.t) if args.s else (s <= args.t), data)
    n_top = float(len(data))
    h = count_freq(data, states)
    print 'Counts'
    print '\t'.join(['state'] + states)
    for s1 in states:
        print '\t'.join([s1] + map(str, [h[s1][s2] for s2 in states]))

    # print '\nGLOBAL STATS', n_total, n_top
    # print '\t'.join(['state'] + states)
    # for s1 in states:
    #     print '\t'.join([s1] + map(str, [global_stats[s1][s2] for s2 in states]))
    print '\nGlobal Enrichment'
    print '\t'.join(['state'] + states)
    for s1 in states:
        print '\t'.join([s1] + map(str, [
            (h[s1][s2] / n_top) / (global_stats[s1][s2] / n_total) if global_stats[s1][s2] > 0 else 0 for s2 in
            states]))
    print '\n per state A->B Enrichment'
    print '\t'.join(['state'] + states)
    for s1 in states:
        n_top_s1 = float(sum(h[s1].values()))
        n_global_s1 = float(sum(global_stats[s1].values()))
        print '\t'.join([s1] + map(str, [
            (h[s1][s2] / n_top_s1) / (global_stats[s1][s2] / n_global_s1) if global_stats[s1][s2] > 0 else 0 for s2 in
            states]))


for fname in args.i:
    compute_state_change_fractions(fname)








