import argparse
import copy
from itertools import izip, combinations
import os
import pprint
import random
import sys
import operator
from utils import open_file, echo, n_combinations, smooth
from metric import *

# from matplotlib import pyplot as plt


def state_transition(probs_a, probs_b):
    a_states = set()
    b_states = set()
    ENRICHMENT_THRESHOLD = 2
    scores = []
    for s in sorted(probs_a, key=state_index):
        score = abs(probs_a[s] - probs_b[s]) * math.log(probs_a[s] / probs_b[s], 2)
        scores.append(score)

        if score > ENRICHMENT_THRESHOLD:
            a_states.add(s)
        elif score < -ENRICHMENT_THRESHOLD:
            b_states.add(s)
    if not a_states:
        a_states.add('A_NONE')
    if not b_states:
        b_states.add('B_NONE')

    return sorted(a_states, key=state_index), sorted(b_states, key=state_index), scores


def get_diff_score(chrom, group_A_segmentations, metric_A, group_B_segmentations, metric_B, BIN_SIZE):
    score = []
    state_transitions = []
    datasets_A = sorted(group_A_segmentations)
    datasets_B = sorted(group_B_segmentations)

    # print chrom
    # iterate over the bins in both segmentation groups
    for a_bins, b_bins in izip(izip(*[iterstate(group_A_segmentations[d][chrom], BIN_SIZE) for d in datasets_A]),
                               izip(*[iterstate(group_B_segmentations[d][chrom], BIN_SIZE) for d in datasets_B])):

        # a_bins contains the states (and bin_starts) from group A
        # b_bins contains the states (and bin_starts) from group B

        # get the posterior probabilities for group A by summing up the
        # the probability vectors for each state
        a_probs = get_prob_vector(a_bins, datasets_A, metric_A)

        # same for group B
        b_probs = get_prob_vector(b_bins, datasets_B, metric_B)

        # compute the score between A and B
        score.append(symmetric_KL_divergence(a_probs, b_probs) ** 2)
        state_transitions.append(state_transition(a_probs, b_probs))
        # print 'A:', a_bins
        # print 'B:', b_bins
        # print

    return score, state_transitions


def iter_posteriors(fname):
    with open_file(fname) as f:
        f.readline()
        f.readline()
        for line in f:
            yield map(float, line.strip().split())


def get_chrom_posteriors_fname(posteriors_dir, fname, chrom):
    basename = os.path.split(fname)[1]
    celltype = basename.split('_segments')[0]
    return os.path.join(posteriors_dir, celltype + '_' + chrom + '_posterior.txt.gz')


def get_diff_score_based_on_posteriors(chrom,
                                       group_A_segmentations,
                                       metric_A,
                                       group_B_segmentations,
                                       metric_B,
                                       BIN_SIZE,
                                       posteriors_dir,
                                       pure_state_transitions):
    score = []
    state_transitions = []
    datasets_A = sorted(group_A_segmentations)
    datasets_B = sorted(group_B_segmentations)

    # print chrom
    # iterate over the bins in both segmentation groups

    a_posteriors_files = [iter_posteriors(get_chrom_posteriors_fname(posteriors_dir, d, chrom)) for d in datasets_A]
    b_posteriors_files = [iter_posteriors(get_chrom_posteriors_fname(posteriors_dir, d, chrom)) for d in datasets_B]

    for (a_bins,
         a_posteriors,
         b_bins,
         b_posteriors) in izip(izip(*[iterstate(group_A_segmentations[d][chrom], BIN_SIZE) for d in datasets_A]),
                               izip(*a_posteriors_files),
                               izip(*[iterstate(group_B_segmentations[d][chrom], BIN_SIZE) for d in datasets_B]),
                               izip(*b_posteriors_files)):

        # a_bins contains the states (and bin_starts) from group A
        # b_bins contains the states (and bin_starts) from group B

        # get the posterior probabilities for group A by summing up the
        # the probability vectors for each state
        a_probs = get_prob_vector_based_on_posteriors(a_posteriors, datasets_A, metric_A)

        # same for group B
        b_probs = get_prob_vector_based_on_posteriors(b_posteriors, datasets_B, metric_B)
        flat_prob_dict = flat_dict({GROUP_A: a_probs, GROUP_B: b_probs})

        closest_state_A, closest_state_B = min([(s1, s2) for s1 in states for s2 in states],
                                                   key=lambda (_s1, _s2):
                                                   KL(flat_prob_dict,
                                                      pure_state_transitions[_s1][_s2][STATE_DISTRIBUTION]))

        pure_state_transitions[closest_state_A][closest_state_B][REAL_COUNTS] += 1

        # compute the score between A and B
        score.append(symmetric_KL_divergence(a_probs, b_probs))
        state_transitions.append(state_transition(a_probs, b_probs))
        # print 'A:', a_bins
        # print 'B:', b_bins
        # print

    for f in a_posteriors_files + b_posteriors_files:
        f.close()

    return score, state_transitions


def iter_state_probs(chrom,
                     group_A_segmentations,
                     metric_A,
                     group_B_segmentations,
                     metric_B,
                     BIN_SIZE,
                     posteriors_dir):

    datasets_A = sorted(group_A_segmentations)
    datasets_B = sorted(group_B_segmentations)

    # print chrom
    # iterate over the bins in both segmentation groups

    a_posteriors_files = [iter_posteriors(get_chrom_posteriors_fname(posteriors_dir, d, chrom)) for d in datasets_A]
    b_posteriors_files = [iter_posteriors(get_chrom_posteriors_fname(posteriors_dir, d, chrom)) for d in datasets_B]

    for (a_bins,
         a_posteriors,
         b_bins,
         b_posteriors) in izip(izip(*[iterstate(group_A_segmentations[d][chrom], BIN_SIZE) for d in datasets_A]),
                               izip(*a_posteriors_files),
                               izip(*[iterstate(group_B_segmentations[d][chrom], BIN_SIZE) for d in datasets_B]),
                               izip(*b_posteriors_files)):

        # a_bins contains the states (and bin_starts) from group A
        # b_bins contains the states (and bin_starts) from group B

        # get the posterior probabilities for group A by summing up the
        # the probability vectors for each state
        a_probs = get_prob_vector_based_on_posteriors(a_posteriors, datasets_A, metric_A)

        # same for group B
        b_probs = get_prob_vector_based_on_posteriors(b_posteriors, datasets_B, metric_B)
        yield {GROUP_A: a_probs,
               GROUP_B: b_probs}

    for f in a_posteriors_files + b_posteriors_files:
        f.close()


def store_output(chrom,
                 signal,
                 state_transitions,
                 out_signal_f,
                 out_states_f,
                 span):
    out_signal_f.write('fixedStep\tchrom=%s\tstart=0\tstep=%d\tspan=%d\n' % (chrom, span, span))
    start = 0
    for bin_no, (bin_signal, (a_states, b_states, scores)) in enumerate(izip(signal, state_transitions)):
        out_signal_f.write('%.2lf\n' % bin_signal)
        out_states_f.write('\t'.join(map(str,
                                         [chrom,
                                          start,
                                          start + BIN_SIZE,
                                          '.',
                                          bin_signal,
                                          '.',
                                          ','.join(a_states) + ' <-> ' + ','.join(b_states),
                                          ','.join('%.1f' % s for s in scores)])) + '\n')

        start += BIN_SIZE


def compute_background_scores(group_A_segmentations,
                              group_B_segmentations,
                              chromosomes,
                              states,
                              posteriors_dir,
                              to_smooth,
                              pure_state_transitions):

    MAX_N_SHUFFLED_SAMPLES = 100
    longest_chromosome = max(chromosomes, key=lambda c: group_A_segmentations.values()[0][c][-1][1])
    print 'Longest chromosome:', longest_chromosome

    n_A_segs = len(group_A_segmentations)
    n_B_segs = len(group_B_segmentations)

    A_segs_keys = tuple(sorted(group_A_segmentations.keys()))
    B_segs_keys = tuple(sorted(group_B_segmentations.keys()))

    # make the union of the two groups
    all_seg_keys = tuple(k for seg_keys in [A_segs_keys, B_segs_keys] for k in seg_keys)
    all_segs = dict(group_A_segmentations.items() + group_B_segmentations.items())

    if n_A_segs == n_B_segs:
        n_comb = n_combinations(len(all_seg_keys), n_A_segs) / 2 - 1
    else:
        n_comb = n_combinations(len(all_seg_keys), n_A_segs) - 1

    def shuffled_samples():

        if n_comb > MAX_N_SHUFFLED_SAMPLES:
            for _ in xrange(MAX_N_SHUFFLED_SAMPLES):
                shuffled_A_segs = tuple(sorted(random.sample(all_seg_keys, n_A_segs)))
                shuffled_B_segs = tuple(sorted(s for s in all_seg_keys if s not in shuffled_A_segs))
                yield shuffled_A_segs, shuffled_B_segs

            pass
        else:
            # generate all shuffled samples
            for shuffled_A_segs in combinations(all_seg_keys, n_A_segs):
                shuffled_B_segs = tuple(s for s in all_seg_keys if s not in shuffled_A_segs)
                yield shuffled_A_segs, shuffled_B_segs

    seen = set([(A_segs_keys, B_segs_keys),
                (B_segs_keys, A_segs_keys)])

    for combo_no, (shuffled_A_segs, shuffled_B_segs) in enumerate(shuffled_samples()):

        if (shuffled_A_segs, shuffled_B_segs) in seen or (shuffled_B_segs, shuffled_A_segs) in seen:
            continue

        seen.add((shuffled_A_segs, shuffled_B_segs))
        seen.add((shuffled_B_segs, shuffled_A_segs))

        print 'Combo no:', combo_no
        print 'Shuffled A:', shuffled_A_segs
        print 'Shuffled B:', shuffled_B_segs
        print

        # reduced_A_segs, reduced_B_segs = pick_bins_at_random(dict((s, all_segs[s]) for s in shuffled_A_segs),
        #                                                      dict((s, all_segs[s]) for s in shuffled_B_segs))

        reduced_A_segs = dict((s, {longest_chromosome: all_segs[s][longest_chromosome]}) for s in shuffled_A_segs)
        reduced_B_segs = dict((s, {longest_chromosome: all_segs[s][longest_chromosome]}) for s in shuffled_B_segs)

        metric_A = learn_metric(reduced_A_segs, states, BIN_SIZE)
        #print_metric(metric_A)

        metric_B = learn_metric(reduced_B_segs, states, BIN_SIZE)
        #print_metric(metric_B)
        for prob_dict in iter_state_probs(longest_chromosome,
                                          reduced_A_segs,
                                          metric_A,
                                          reduced_B_segs,
                                          metric_B,
                                          BIN_SIZE,
                                          posteriors_dir):

            flat_prob_dict = flat_dict(prob_dict)
            closest_state_A, closest_state_B = min([(s1, s2) for s1 in states for s2 in states],
                                                   key=lambda (_s1, _s2):
                                                   KL(flat_prob_dict,
                                                      pure_state_transitions[_s1][_s2][STATE_DISTRIBUTION]))

            pure_state_transitions[closest_state_A][closest_state_B][SHUFFLED_FREQUENCY] += 1
            pure_state_transitions[closest_state_B][closest_state_A][SHUFFLED_FREQUENCY] += 1


    # normalize the shuffled frequencies
    total = sum(pure_state_transitions[s1][s2][SHUFFLED_FREQUENCY] for s1 in states for s2 in states)
    for s1 in states:
        for s2 in states:
            pure_state_transitions[s1][s2][SHUFFLED_FREQUENCY] /= float(total)


parser = argparse.ArgumentParser()

parser.add_argument('-a', nargs='+', help='segmentation files of group A')
parser.add_argument('-b', nargs='+', help='segmentation files of group B')
parser.add_argument('-p', '--posteriors', help='directory with posteriors')

parser.add_argument('--bin-size', type=int, help='Bin size, default: %(default)s', default=200)
parser.add_argument('-f', '--fdr-threshold', type=float, help='FDR threshold: %(default)s', default=0.01)

parser.add_argument('-s', '--smooth', action='store_true', help='smooth signal')
parser.add_argument('-o', help='output prefix')


args = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    exit(1)

BIN_SIZE = args.bin_size
posteriors_dir = args.posteriors
echo('Smoothing:', args.smooth)
echo('Posteriors dir:', posteriors_dir)

group_A_segmentations, states = read_segmentations(args.a)
group_B_segmentations, _ = read_segmentations(args.b)

chromosomes = sorted(reduce(operator.and_, [set(s.keys())
                                         for segs in [group_A_segmentations, group_B_segmentations]
                                         for s in segs.values()]))
# print chromosomes
group_A_segmentations = filter_chroms(group_A_segmentations, chromosomes)
group_B_segmentations = filter_chroms(group_B_segmentations, chromosomes)

metric_A = learn_metric(group_A_segmentations, states, BIN_SIZE)
print_metric(metric_A)

metric_B = learn_metric(group_B_segmentations, states, BIN_SIZE)
print_metric(metric_B)
print_average_metric(states, metric_A, metric_B)

fdr_threshold = args.fdr_threshold

pure_state_transitions = generate_pure_state_transitions(states, metric_A, metric_B)

if fdr_threshold < 1:

    echo('Learning significance threshold for at p-value:', fdr_threshold)
    compute_background_scores(group_A_segmentations,
                              group_B_segmentations,
                              chromosomes,
                              states,
                              posteriors_dir,
                              args.smooth,
                              pure_state_transitions)

# longest_chromosome = max(chromosomes, key=lambda c: group_A_segmentations.values()[0][c][-1][1])
# print 'significance threshold at', fdr_threshold, ':', score_cutoff
# exit(1)

out_fname = args.o


# threshold = determine_threshold(group_A_segmentations, metric_A,
#                                 group_B_segmentations, metric_B,
#                                 fdr_threshold=.01)
real_scores = {}
echo('Storing output:', out_fname)
total_genomic_bins = 0
with open_file(out_fname + '.wig.gz', 'w') as out_signal_f, \
     open_file(out_fname + '.bed.gz', 'w') as out_states_f:
    title = os.path.split(out_fname)[1].replace('.wig', '').replace('.gz', '')
    out_signal_f.write('track type=wiggle_0 name="%s" description="%s"\n' % (title, title))
    for chrom in chromosomes:
        echo('Chromosome:', chrom)
        signal, state_transitions = get_diff_score_based_on_posteriors(chrom,
                                                                       group_A_segmentations,
                                                                       metric_A,
                                                                       group_B_segmentations,
                                                                       metric_B,
                                                                       BIN_SIZE,
                                                                       posteriors_dir,
                                                                       pure_state_transitions)
        if args.smooth:
            signal = smooth(signal)

        real_scores[chrom] = signal
        total_genomic_bins += len(signal)
        # if chrom == longest_chromosome:
        #
        #     hist_info = plt.hist([bgr_scores, signal], bins=50, color=['blue', 'red'], label=['bgr', 'real'],
        #                          normed=1)
        #     max_y = max(v for a in [hist_info[0][0][2:], hist_info[0][1][2:]] for v in a)

        store_output(chrom,
                     signal,
                     state_transitions,
                     out_signal_f=out_signal_f,
                     out_states_f=out_states_f,
                     span=BIN_SIZE)

if fdr_threshold >= 1:
    exit()


# for s1 in states:
#     for s2 in states:
#         p = pure_state_transitions[s1][s2][SHUFFLED_FREQUENCY]
#         pure_state_transitions[s1][s2][P_VALUE] = 1 - normal_cdf(pure_state_transitions[s1][s2][REAL_COUNTS],
#                                                                  total_genomic_bins * p,
#                                                                  math.sqrt(total_genomic_bins * p * (1 - p)))

print_pure_state_transitions(pure_state_transitions)

#
# hist_info = plt.hist([background_scores, all_scores], bins=50, color=['blue', 'red'], label=['bgr', 'real'],
#                      normed=1)
# max_y = max(v for a in [hist_info[0][0][2:], hist_info[0][1][2:]] for v in a)
#
# xmin, xmax = plt.xlim()
# plt.xlim(1, xmax)
# plt.ylim(0, max_y)
# plt.legend()
# plt.savefig(out_fname + '_score_histogram.png')
#
#
