import argparse, sys
from utils import *
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

def read_scores(fname):
    echo("Reading:", fname)

    scores = []
    with open_file(fname) as in_f:
        for l in in_f:
            buf = l.split()
            key = '\t'.join(buf[:3])
            scores.append((abs(float(buf[4])), key))

    return sorted(scores, reverse=True)



parser = argparse.ArgumentParser()

parser.add_argument('-i', nargs='+', help='scores')
parser.add_argument('-o', help='output png')


args = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    exit(1)

scores = {}
for fname in args.i:
    # mean_distance_matrix.split_0_A.overall_score.summary.bed.gz
    parts = fname.split('.split_')
    method_name = parts[0]
    split_no, split_part = parts[1].split('.')[0].split('_')

    print fname
    print method_name, split_no, split_part

    if method_name not in scores:
        scores[method_name] = {}

    if split_no not in scores[method_name]:
        scores[method_name][split_no] = []

    scores[method_name][split_no].append(fname)

echo('Computing rank overlaps')

plt.figure(figsize=(20, 10))

N_POINTS = 150000

MDM = 'mean_distance_matrix'
PRDM = 'per_replicate_distance_matrix'


mean_scores = {MDM: [0] * N_POINTS,
               PRDM: [0] * N_POINTS }

mean_squared = {MDM: [0] * N_POINTS,
                PRDM: [0] * N_POINTS }


n_samples = {MDM: 0, PRDM: 0}
for method_name in scores:
    for split_no in scores[method_name]:
        echo('Plotting', method_name, split_no)
        n_samples[method_name] += 1

        if len(scores[method_name][split_no]) != 2:
            print 'Error:', method_name, split_no
            exit(1)

        a_fname, b_fname = scores[method_name][split_no]

        a_scores = read_scores(a_fname)
        b_scores = read_scores(b_fname)


        a_scores = a_scores[:N_POINTS]
        b_scores = b_scores[:N_POINTS]

        # print a_scores[:5]
        # print b_scores[:5]

        n_scores = len(a_scores)

        overlap_fraction = [0] * n_scores

        c_overlap = 0

        a_set = set()
        b_set = set()

        for rank, ((s1, key_1), (s2, key_2)) in enumerate(izip(a_scores, b_scores)):

            a_set.add(key_1)
            b_set.add(key_2)

            if key_1 == key_2:
                c_overlap += 1
            else:
                if key_2 in a_set:
                    c_overlap += 1

                if key_1 in b_set:
                    c_overlap += 1

            overlap_fraction[rank] = float(c_overlap) / (rank + 1)

            mean_scores[method_name][rank] += overlap_fraction[rank]
            mean_squared[method_name][rank] += overlap_fraction[rank] ** 2

        plt.plot(overlap_fraction, color='blue' if method_name == 'mean_distance_matrix' else 'red') # , label=(method_name + ' ' + split_no))

print n_samples
# plt.legend(loc='lower right')
plt.savefig(args.o + '.individual.png')
import numpy as np

plt.figure(figsize=(20, 10))
for method_name in mean_scores:
    color = 'blue' if method_name == 'mean_distance_matrix' else 'red'
    means = np.array([v / n_samples[method_name] for v in mean_scores[method_name]])
    stds = np.array([math.sqrt(v / n_samples[method_name] - m ** 2) for v, m in izip(mean_squared[method_name], means)])

    plt.plot(range(N_POINTS), means, color=color, label=method_name)
    plt.fill_between(range(N_POINTS), means - stds, means + stds, color=color, alpha=0.2)

plt.legend()
plt.savefig(args.o + '.aggregate.png')


