import argparse
import gzip

import sys
import math

def open_file(fname, mode='r'):
    return gzip.open(fname, mode) if fname.endswith('.gz') else open(fname, mode)

def iter_posteriors(fname):
    with open_file(fname) as f:
        f.readline()
        f.readline()
        for line in f:
            yield tuple(map(lambda v: round(float(v), 3), line.strip().split()))


def entropy(d):
    e = 0
    for p in d.values():
        if p > 0:
            e -= p * math.log(p, 2)
    return e


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-a', nargs='+', help='posterior files')

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)

    posterior_squared = None
    states = None

    for fname in args.a:
        with open_file(fname) as f:
            print >>sys.stderr, fname
            f.readline()
            states = f.readline().strip().split('\t')
            if posterior_squared is None:
                posterior_squared = dict((s1, dict((s2, 0) for s2 in states)) for s1 in states)

            for line in f:
                posteriors = map(lambda v: round(float(v), 3), line.strip().split())
                non_zero = [(s, p) for s, p in zip(states, posteriors) if p > 0]
                for s1, p1 in non_zero:
                    for s2, p2 in non_zero:
                        posterior_squared[s1][s2] += p1 * p2
                # for i, s1 in enumerate(states):
                #     if posteriors[i] > 0:
                #         for j, s2 in enumerate(states):
                #             posterior_squared[s1][s2] += posteriors[i] * posteriors[j]
    print '\t'.join(['state'] + states + ['Entropy'])
    for s1 in states:
        total = sum(posterior_squared[s1].values())
        print s1,
        for s2 in states:
            posterior_squared[s1][s2] /= total
            print '\t' + str(round(100 * posterior_squared[s1][s2], 2)),

        print '\t' + '%.2lf' % (entropy(posterior_squared[s1]))


