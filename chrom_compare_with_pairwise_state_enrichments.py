import argparse
from itertools import izip, combinations
import os
import random
import sys
import operator
import itertools
from utils import open_file, echo, n_combinations, smooth
from metric import *


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


def get_diff_score(chrom,
                   group_A_segmentations,
                   metric_A,
                   group_B_segmentations,
                   metric_B,
                   BIN_SIZE,
                   pure_state_transitions,
                   closest_transition_cache,
                   scores_cache):

    score = []
    state_transitions = []

    datasets_A = sorted(group_A_segmentations)
    datasets_B = sorted(group_B_segmentations)

    # print chrom
    state_pairs = [(s1, s2) for s1 in states for s2 in states]

    # iterate over the bins in both segmentation groups
    for a_bins, b_bins in izip(izip(*[iterstate(group_A_segmentations[d][chrom], BIN_SIZE) for d in datasets_A]),
                               izip(*[iterstate(group_B_segmentations[d][chrom], BIN_SIZE) for d in datasets_B])):

        # a_bins contains the states (and bin_starts) from group A
        # b_bins contains the states (and bin_starts) from group B

        # get the posterior probabilities for group A by summing up the
        # the probability vectors for each state
        a_probs = get_prob_vector(a_bins, datasets_A, metric_A, states)

        # same for group B
        b_probs = get_prob_vector(b_bins, datasets_B, metric_B, states)

        key = (tuple([_s for _, _s in a_bins]),
               tuple([_s for _, _s in b_bins]))

        if key not in closest_transition_cache:
            flat_prob_dict = dict((group + ' ' + key, d[key])
                                  for group, d in [('A', a_probs), ('B', b_probs)] for key in d)

            closest_transition_cache[key] = min(state_pairs,
                                                key=lambda (_s1, _s2):
                                                KL(flat_prob_dict,
                                                   pure_state_transitions[_s1][_s2][STATE_DISTRIBUTION]))

            scores_cache[key] = symmetric_KL_divergence(a_probs, b_probs)

        closest_state_A, closest_state_B = closest_transition_cache[key]

        pure_state_transitions[closest_state_A][closest_state_B][REAL_COUNTS] += 1

        # compute the score between A and B
        score.append(scores_cache[key])
        state_transitions.append((closest_state_A,
                                  closest_state_B,
                                  [s for _, s in a_bins],
                                  [s for _, s in b_bins]))
        # print 'A:', a_bins
        # print 'B:', b_bins
        # print

    return score, state_transitions

def store_output(chrom,
                 signal,
                 state_transitions,
                 out_signal_f,
                 out_states_f,
                 span):
    out_signal_f.write('fixedStep\tchrom=%s\tstart=0\tstep=%d\tspan=%d\n' % (chrom, span, span))
    start = 0
    for bin_no, (bin_signal, (closest_state_A,
                              closest_state_B,
                              a_states,
                              b_states)) in enumerate(izip(signal, state_transitions)):

        out_signal_f.write('%.2lf\n' % bin_signal)
        out_states_f.write('\t'.join(map(str,
                                         [chrom,
                                          start,
                                          start + BIN_SIZE,
                                          '.',
                                          bin_signal,
                                          '.',
                                          closest_state_A + '-' + closest_state_B,
                                          ','.join(a_states) + ' <-> ' + ','.join(b_states)
                                          ])) + '\n')

        start += BIN_SIZE


def compute_background_frequencies(group_A_segmentations,
                                   group_B_segmentations,
                                   chromosomes,
                                   states):

    MAX_N_SHUFFLED_SAMPLES = 100
    # longest_chromosome = max(chromosomes, key=lambda c: group_A_segmentations.values()[0][c][-1][1])
    # print 'Longest chromosome:', longest_chromosome

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

    combo_no = 0

    shuffled_frequencies = dict((s1, dict((s2, 0) for s2 in states)) for s1 in states)
    background_scores = []

    state_pairs = [(s1, s2) for s1 in states for s2 in states]

    for shuffled_A_datasets, shuffled_B_datasets in shuffled_samples():

        if (shuffled_A_datasets, shuffled_B_datasets) in seen or (shuffled_B_datasets, shuffled_A_datasets) in seen:
            continue

        seen.add((shuffled_A_datasets, shuffled_B_datasets))
        seen.add((shuffled_B_datasets, shuffled_A_datasets))

        combo_no += 1
        echo('Combo no:', combo_no)
        print 'Shuffled A:', shuffled_A_datasets
        print 'Shuffled B:', shuffled_B_datasets
        print

        # reduced_A_segs, reduced_B_segs = pick_bins_at_random(dict((s, all_segs[s]) for s in shuffled_A_segs),
        #                                                      dict((s, all_segs[s]) for s in shuffled_B_segs))

        # reduced_A_segs = dict((s, {longest_chromosome: all_segs[s][longest_chromosome]}) for s in shuffled_A_segs)
        # reduced_B_segs = dict((s, {longest_chromosome: all_segs[s][longest_chromosome]}) for s in shuffled_B_segs)

        shuffled_A_segmentations = dict((dataset, all_segs[dataset]) for dataset in shuffled_A_datasets)
        shuffled_B_segmentations = dict((dataset, all_segs[dataset]) for dataset in shuffled_B_datasets)

        metric_A = learn_metric(shuffled_A_segmentations, states, BIN_SIZE)
        #print_metric(metric_A)

        metric_B = learn_metric(shuffled_B_segmentations, states, BIN_SIZE)
        #print_metric(metric_B)
        # print '*' * 50
        # print_average_metric(states, metric_A, metric_B)
        # print '*' * 50
        current_pure_state_transitions = generate_pure_state_transitions(states, metric_A, metric_B)

        # add the pseudo-counts
        for s1 in shuffled_frequencies:
            for s2 in shuffled_frequencies[s1]:
                shuffled_frequencies[s1][s2] += PSEUDO_COUNT

        closest_transition_cache = {}
        scores_cache = {}

        datasets_A = sorted(shuffled_A_segmentations)
        datasets_B = sorted(shuffled_B_segmentations)

        # print chrom
        # iterate over the bins in both segmentation groups
        for chrom in chromosomes:
            # echo('Chromosome:', chrom)
            for (a_bins, b_bins) in izip(izip(*[iterstate(shuffled_A_segmentations[d][chrom], BIN_SIZE) for d in datasets_A]),
                                         izip(*[iterstate(shuffled_B_segmentations[d][chrom], BIN_SIZE) for d in datasets_B])):

                # a_bins contains the states (and bin_starts) from group A
                # b_bins contains the states (and bin_starts) from group B

                key = (tuple([_s for _, _s in a_bins]),
                       tuple([_s for _, _s in b_bins]))

                if key not in closest_transition_cache:
                    a_probs = get_prob_vector(a_bins, datasets_A, metric_A, states)

                    # same for group B
                    b_probs = get_prob_vector(b_bins, datasets_B, metric_B, states)
                    flat_prob_dict = dict((group + ' ' + key, d[key])
                                          for group, d in [('A', a_probs), ('B', b_probs)] for key in d)

                    closest_transition_cache[key] = min(state_pairs,
                                                        key=lambda (_s1, _s2):
                                                        KL(flat_prob_dict,
                                                           current_pure_state_transitions[_s1][_s2][STATE_DISTRIBUTION]))

                    scores_cache[key] = symmetric_KL_divergence(a_probs, b_probs)

                closest_state_A, closest_state_B = closest_transition_cache[key]
                background_scores.append(scores_cache[key])
                # combos_to_check = list(itertools.product([state for _, state in a_states_set],
                #                                          [state for _, state in b_states_set]))

                # print a_states_set, b_states_set
                # print list(combos_to_check)
                # exit(1)
                # get the posterior probabilities for group A by summing up the
                # the probability vectors for each state

                shuffled_frequencies[closest_state_A][closest_state_B] += 1
                shuffled_frequencies[closest_state_B][closest_state_A] += 1

    # normalize the shuffled frequencies
    #     print _n, _t
    # exit(1)
    total = sum(shuffled_frequencies[s1][s2] for s1 in states for s2 in states)
    for s1 in states:
        for s2 in states:
            shuffled_frequencies[s1][s2] /= float(total)

    return shuffled_frequencies, sorted(background_scores)

parser = argparse.ArgumentParser()

parser.add_argument('-a', nargs='+', help='segmentation files of group A')
parser.add_argument('-b', nargs='+', help='segmentation files of group B')
parser.add_argument('--bin-size', type=int, help='Bin size, default: %(default)s', default=200)
parser.add_argument('-f', '--fdr-threshold', type=float, help='FDR threshold: %(default)s', default=0.01)

parser.add_argument('-o', help='output prefix')
parser.add_argument('-s', '--smooth', action='store_true', help='smooth signal')


args = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    exit(1)

BIN_SIZE = args.bin_size

echo('Smoothing:', args.smooth)

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

echo('Learning significance threshold for at p-value:', fdr_threshold)
shuffled_frequencies, background_scores = compute_background_frequencies(group_A_segmentations,
                                                                         group_B_segmentations,
                                                                         chromosomes,
                                                                         states)
print SHUFFLED_FREQUENCY
print '\t'.join(['State'] + states)
for s1 in sorted(shuffled_frequencies, key=state_key):
    print '\t'.join([s1] + map(str, [shuffled_frequencies[s1][s2] for s2 in states]))
print

for s1 in shuffled_frequencies:
    for s2 in shuffled_frequencies[s1]:
        pure_state_transitions[s1][s2][SHUFFLED_FREQUENCY] = shuffled_frequencies[s1][s2]

    # longest_chromosome = max(chromosomes, key=lambda c: group_A_segmentations.values()[0][c][-1][1])
    # print 'significance threshold at', fdr_threshold, ':', score_cutoff
    # exit(1)

out_fname = args.o

# threshold = determine_threshold(group_A_segmentations, metric_A,
#                                 group_B_segmentations, metric_B,
#                                 fdr_threshold=.01)
echo('Storing output')
real_scores = {}
total_genomic_bins = 0
with open_file(out_fname + '.wig.gz', 'w') as out_signal_f, \
     open_file(out_fname + '.bed.gz', 'w') as out_states_f:

    title = os.path.split(out_fname)[1].replace('.wig', '').replace('.gz', '')
    out_signal_f.write('track type=wiggle_0 name="%s" description="%s"\n' % (title, title))
    closest_transition_cache = {}
    scores_cache = {}
    for chrom in chromosomes:
        echo('Chromosome:', chrom)
        signal, state_transitions = get_diff_score(chrom,
                                                   group_A_segmentations,
                                                   metric_A,
                                                   group_B_segmentations,
                                                   metric_B,
                                                   BIN_SIZE,
                                                   pure_state_transitions,
                                                   closest_transition_cache,
                                                   scores_cache)
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

#
# for s1 in states:
#     for s2 in states:
#         p = pure_state_transitions[s1][s2][SHUFFLED_FREQUENCY]
#         pure_state_transitions[s1][s2][P_VALUE] = 1 - normal_cdf(pure_state_transitions[s1][s2][REAL_COUNTS],
#                                                                  total_genomic_bins * p,
#                                                                  math.sqrt(total_genomic_bins * p * (1 - p)))

for s1 in states:
    for s2 in states:
        pure_state_transitions[s1][s2][REAL_FREQUENCY] = float(pure_state_transitions[s1][s2][REAL_COUNTS]) / total_genomic_bins
        pure_state_transitions[s1][s2][P_VALUE] = binom_test(pure_state_transitions[s1][s2][REAL_COUNTS],
                                                             total_genomic_bins,
                                                             expected_freq=shuffled_frequencies[s1][s2])
print_pure_state_transitions(pure_state_transitions)
print 'Number of state transitions:', len(states) ** 2
echo('Done')