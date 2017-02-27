import argparse
from itertools import izip, combinations
import os
import random
import sys
import operator
import pickle
from permute_marks_in_binarized_ChromHMM import compute_background_scores_by_shuffling_marks
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

OVERALL_SCORE = 'overall_score'
def get_overall_and_per_state_diff_score(chrom,
                                           group_A_segmentations,
                                           metric_A,
                                           group_B_segmentations,
                                           metric_B,
                                           BIN_SIZE,
                                           states,
                                           scores_cache,
                                           pure_state_transitions,
                                           closest_transition_cache):

    scores = dict((k, []) for k in [OVERALL_SCORE] + states)

    state_transitions = []

    state_pairs = [(s1, s2) for s1 in states for s2 in states]

    datasets_A = sorted(group_A_segmentations)
    datasets_B = sorted(group_B_segmentations)

    # iterate over the bins in both segmentation groups
    for a_bins, b_bins in izip(izip(*[iterstate(group_A_segmentations[d][chrom], BIN_SIZE) for d in datasets_A]),
                               izip(*[iterstate(group_B_segmentations[d][chrom], BIN_SIZE) for d in datasets_B])):

        # get the posterior probabilities for group A by summing up the
        # the probability vectors for each state
        a_probs = get_prob_vector(a_bins, datasets_A, metric_A, states)

        # same for group B
        b_probs = get_prob_vector(b_bins, datasets_B, metric_B, states)

        key = (tuple([_s for _, _s in a_bins]),
               tuple([_s for _, _s in b_bins]))

        if key not in scores_cache:
            scores_cache[key] = dict((k, 0) for k in [OVERALL_SCORE] + states)

            scores_cache[key][OVERALL_SCORE] = symmetric_KL_divergence(a_probs, b_probs)

            for state in states:
                scores_cache[key][state] = a_probs[state] - b_probs[state]

            flat_prob_dict = dict((group + ' ' + key, d[key])
                                  for group, d in [('A', a_probs), ('B', b_probs)] for key in d)

            closest_transition_cache[key] = min(state_pairs,
                                                key=lambda (_s1, _s2):
                                                KL(flat_prob_dict,
                                                   pure_state_transitions[_s1][_s2][STATE_DISTRIBUTION]))

        # compute the score between A and B
        for score_type in scores_cache[key]:
            scores[score_type].append(scores_cache[key][score_type])

        closest_state_A, closest_state_B = closest_transition_cache[key]

        pure_state_transitions[closest_state_A][closest_state_B][REAL_COUNTS] += 1

        # compute the score between A and B
        state_transitions.append((closest_state_A,
                                  closest_state_B,
                                  [s for _, s in a_bins],
                                  [s for _, s in b_bins]))

    return scores, state_transitions


def store_wig(chrom,
             signal,
             out_signal_f,
             span):

    out_signal_f.write('fixedStep\tchrom=%s\tstart=0\tstep=%d\tspan=%d\n' % (chrom, span, span))

    for bin_no, bin_signal in enumerate(signal):

        out_signal_f.write('%.2lf\n' % bin_signal)




def store_bed(chrom,
              state_transitions,
              signal,
              score_to_fdr,
              out_states_f):

    start = 0
    for bin_no, (bin_signal, (closest_state_A,
                              closest_state_B,
                              a_states,
                              b_states)) in enumerate(izip(signal, state_transitions)):

        out_states_f.write('\t'.join(map(str,
                                         [chrom,
                                          start,
                                          start + BIN_SIZE,
                                          '.',
                                          bin_signal,
                                          '.',
                                          closest_state_A + '-' + closest_state_B,
                                          ','.join(a_states) + '/' + ','.join(b_states),
                                          score_to_fdr[bin_signal]
                                          ])) + '\n')

        start += BIN_SIZE

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
                                          ','.join(a_states) + '/' + ','.join(b_states)
                                          ])) + '\n')

        start += BIN_SIZE


def compute_background_scores_by_shuffling_segmentations(group_A_segmentations,
                                                         group_B_segmentations,
                                                         chromosomes,
                                                         states,
                                                         metric_A,
                                                         metric_B,
                                                         n_perm=100,
                                                         to_smooth=False
                                                         ):

    # longest_chromosome = max(chromosomes, key=lambda c: group_A_segmentations.values()[0][c][-1][1])
    # print 'Longest chromosome:', longest_chromosome

    n_A_segs = len(group_A_segmentations)
    n_B_segs = len(group_B_segmentations)

    A_segs_keys = tuple(sorted(group_A_segmentations.keys()))
    B_segs_keys = tuple(sorted(group_B_segmentations.keys()))

    # make the union of the two groups
    all_seg_keys = tuple(k for seg_keys in [A_segs_keys, B_segs_keys] for k in seg_keys)
    all_segs = dict(group_A_segmentations.items() + group_B_segmentations.items())
    all_metrics = dict(metric_A.items() + metric_B.items())

    if n_A_segs == n_B_segs:
        n_comb = n_combinations(len(all_seg_keys), n_A_segs) / 2 - 1
    else:
        n_comb = n_combinations(len(all_seg_keys), n_A_segs) - 1

    def shuffled_samples(n_perm):

        if args.all_for_null:
            seen = set()
        else:
            seen = set([(A_segs_keys, B_segs_keys),
                        (B_segs_keys, A_segs_keys)])

        if n_comb > n_perm:

            if args.all_for_null:
                yield A_segs_keys, B_segs_keys
                seen = set([(A_segs_keys, B_segs_keys),
                            (B_segs_keys, A_segs_keys)])
                n_perm -= 1

            for _ in xrange(n_perm):

                shuffled_A_segs = tuple(sorted(random.sample(all_seg_keys, n_A_segs)))
                shuffled_B_segs = tuple(sorted(s for s in all_seg_keys if s not in shuffled_A_segs))

                while (shuffled_A_segs, shuffled_B_segs) in seen and (shuffled_B_segs, shuffled_A_segs) in seen:
                    shuffled_A_segs = tuple(sorted(random.sample(all_seg_keys, n_A_segs)))
                    shuffled_B_segs = tuple(sorted(s for s in all_seg_keys if s not in shuffled_A_segs))

                seen.add((shuffled_A_segs, shuffled_B_segs))
                seen.add((shuffled_B_segs, shuffled_A_segs))

                yield shuffled_A_segs, shuffled_B_segs

        else:
            # generate all shuffled samples
            for shuffled_A_segs in combinations(all_seg_keys, n_A_segs):
                shuffled_B_segs = tuple(s for s in all_seg_keys if s not in shuffled_A_segs)

                if (shuffled_A_segs, shuffled_B_segs) in seen or (shuffled_B_segs, shuffled_A_segs) in seen:
                    continue

                seen.add((shuffled_A_segs, shuffled_B_segs))
                seen.add((shuffled_B_segs, shuffled_A_segs))

                yield shuffled_A_segs, shuffled_B_segs

    combo_no = 0

    # shuffled_frequencies = dict((s1, dict((s2, 0) for s2 in states)) for s1 in states)
    background_scores = []

    # state_pairs = [(s1, s2) for s1 in states for s2 in states]

    for shuffled_A_datasets, shuffled_B_datasets in shuffled_samples(n_perm):

        combo_no += 1
        echo('Combo no:', combo_no)
        print 'Shuffled A:', shuffled_A_datasets
        print 'Shuffled B:', shuffled_B_datasets
        print

        # reduced_A_segs, reduced_B_segs = pick_bins_at_random(dict((s, all_segs[s]) for s in shuffled_A_segs),
        #                                                      dict((s, all_segs[s]) for s in shuffled_B_segs))

        # reduced_A_segs = dict((s, {longest_chromosome: all_segs[s][longest_chromosome]}) for s in shuffled_A_segs)
        # reduced_B_segs = dict((s, {longest_chromosome: all_segs[s][longest_chromosome]}) for s in shuffled_B_segs)

        # N_RANDOM_LOCATIONS = 100
        # N_LOCUS_LENGTH = 1000000
        # shuffled_A_segmentations, shuffled_B_segmentations = sample_locations(dict((dataset, all_segs[dataset]) for dataset in shuffled_A_datasets),
        #                                                                       dict((dataset, all_segs[dataset]) for dataset in shuffled_B_datasets),
        #                                                                       N_RANDOM_LOCATIONS,
        #                                                                       N_LOCUS_LENGTH)

        shuffled_A_segmentations = dict((dataset, all_segs[dataset]) for dataset in shuffled_A_datasets)
        shuffled_B_segmentations = dict((dataset, all_segs[dataset]) for dataset in shuffled_B_datasets)

        shuffled_metric_A = dict((dataset, all_metrics[dataset]) for dataset in shuffled_A_datasets)
        shuffled_metric_B = dict((dataset, all_metrics[dataset]) for dataset in shuffled_B_datasets)

        # shuffled_metric_A = learn_metric(shuffled_A_segmentations, states, BIN_SIZE)
        # print_metric(shuffled_metric_A)
        #
        # shuffled_metric_B = learn_metric(shuffled_B_segmentations, states, BIN_SIZE)

        # print_metric(shuffled_metric_B)
        # print '*' * 50
        # print_average_metric(states, shuffled_metric_A, shuffled_metric_B)
        # print '*' * 50
        current_pure_state_transitions = generate_pure_state_transitions(states, shuffled_metric_A, shuffled_metric_B)

        # add the pseudo-counts
        # for s1 in shuffled_frequencies:
        #     for s2 in shuffled_frequencies[s1]:
        #         shuffled_frequencies[s1][s2] += PSEUDO_COUNT

        # closest_transition_cache = {}
        scores_cache = {}

        datasets_A = sorted(shuffled_A_segmentations)
        datasets_B = sorted(shuffled_B_segmentations)

        # print chrom
        # iterate over the bins in both segmentation groups
        # for chrom in random.sample(chromosomes, 2):
        for chrom in chromosomes:
            echo('Chromosome:', chrom)
            chrom_background_scores = []
            for (a_bins, b_bins) in izip(izip(*[iterstate(shuffled_A_segmentations[d][chrom], BIN_SIZE) for d in datasets_A]),
                                         izip(*[iterstate(shuffled_B_segmentations[d][chrom], BIN_SIZE) for d in datasets_B])):

                # a_bins contains the states (and bin_starts) from group A
                # b_bins contains the states (and bin_starts) from group B

                key = (tuple([_s for _, _s in a_bins]),
                       tuple([_s for _, _s in b_bins]))

                if key not in scores_cache:
                    a_probs = get_prob_vector(a_bins, datasets_A, shuffled_metric_A, states)

                    # same for group B
                    b_probs = get_prob_vector(b_bins, datasets_B, shuffled_metric_B, states)
                    # flat_prob_dict = dict((group + ' ' + key, d[key])
                    #                       for group, d in [('A', a_probs), ('B', b_probs)] for key in d)

                    # closest_transition_cache[key] = min(state_pairs,
                    #                                     key=lambda (_s1, _s2):
                    #                                     KL(flat_prob_dict,
                    #                                        current_pure_state_transitions[_s1][_s2][STATE_DISTRIBUTION]))

                    scores_cache[key] = symmetric_KL_divergence(a_probs, b_probs)

                # closest_state_A, closest_state_B = closest_transition_cache[key]
                chrom_background_scores.append(scores_cache[key])
                # combos_to_check = list(itertools.product([state for _, state in a_states_set],
                #                                          [state for _, state in b_states_set]))

                # print a_states_set, b_states_set
                # print list(combos_to_check)
                # exit(1)
                # get the posterior probabilities for group A by summing up the
                # the probability vectors for each state

                # shuffled_frequencies[closest_state_A][closest_state_B] += 1
                # shuffled_frequencies[closest_state_B][closest_state_A] += 1

            if to_smooth:
                chrom_background_scores = smooth(chrom_background_scores)

            background_scores.extend(chrom_background_scores)

    # normalize the shuffled frequencies
    #     print _n, _t
    # exit(1)
    # total = sum(shuffled_frequencies[s1][s2] for s1 in states for s2 in states)
    # for s1 in states:
    #     for s2 in states:
    #         shuffled_frequencies[s1][s2] /= float(total)

    return sorted(background_scores)


parser = argparse.ArgumentParser()

parser.add_argument('-a', nargs='+', help='segmentation files of group A')
parser.add_argument('-b', nargs='+', help='segmentation files of group B')
parser.add_argument('--bin-size', type=int, help='Bin size, default: %(default)s', default=200)
parser.add_argument('-f', '--fdr-threshold', type=float, help='FDR threshold: %(default)s', default=0.01)

parser.add_argument('-o', '--output', help='output prefix')
parser.add_argument('--background-scores', help='bacground_scores.pickle')

parser.add_argument('--chrom-hmm-binarized', nargs='+', help='path to ChromHMM binarized files to compute background scores')
parser.add_argument('--chrom-hmm-model-path', help='path to the ChromHMM model to use for shuffled marks')

parser.add_argument('-s', '--smooth', action='store_true', help='smooth signal')
parser.add_argument('--use-closest-rep', action='store_true', help='use only closest replicate to learn the metric')
parser.add_argument('--all-for-null', action='store_true', default=False, help='include the original combo in the background model')


args = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    exit(1)

BIN_SIZE = args.bin_size
echo('cmd:', ' '.join(sys.argv))
for arg in sorted(vars(args)):
    echo(arg , '=', getattr(args, arg))

echo('Smoothing:', args.smooth)
use_all_reps = not args.use_closest_rep
if use_all_reps:
    learn_metric = learn_metric_from_all_replicates

group_A_segmentations, states = read_segmentations(args.a)
group_B_segmentations, _ = read_segmentations(args.b)

chromosomes = sorted(reduce(operator.and_, [set(s.keys())
                                         for segs in [group_A_segmentations, group_B_segmentations]
                                         for s in segs.values()]),
                     key=lambda c: int(c.replace('chr','')) if re.search(r'^chr(\d+)$', c) else 100)
# print chromosomes
group_A_segmentations = filter_chroms(group_A_segmentations, chromosomes)
group_B_segmentations = filter_chroms(group_B_segmentations, chromosomes)

metric_A = learn_metric(group_A_segmentations, states, BIN_SIZE, other_group=group_B_segmentations)
print_metric(metric_A)

metric_B = learn_metric(group_B_segmentations, states, BIN_SIZE, other_group=group_A_segmentations)
print_metric(metric_B)
print_average_metric(states, metric_A, metric_B)

fdr_threshold = args.fdr_threshold

if fdr_threshold < 1:
    MAX_N_SHUFFLED_SAMPLES = 100
    if args.background_scores:
        echo('Loading background scores from :', args.background_scores)
        background_scores = pickle.load(open(args.background_scores))
    elif args.chrom_hmm_binarized:
        background_scores = compute_background_scores_by_shuffling_marks(args.chrom_hmm_binarized,
                                                                         args.chrom_hmm_model_path,
                                                                         n_perms=MAX_N_SHUFFLED_SAMPLES,
                                                                         n_group_A=len(group_A_segmentations),
                                                                         BIN_SIZE=BIN_SIZE,
                                                                         to_smooth=args.smooth)

    else:
        echo('Learning significance threshold for at p-value:', fdr_threshold)
        background_scores = compute_background_scores_by_shuffling_segmentations(group_A_segmentations,
                                                                                 group_B_segmentations,
                                                                                 chromosomes[:1], #['chr1'],
                                                                                 states,
                                                                                 metric_A,
                                                                                 metric_B,
                                                                                 n_perm=MAX_N_SHUFFLED_SAMPLES,
                                                                                 to_smooth=args.smooth)


out_fname = args.output

real_scores = {}
real_state_transitions = {}
echo('Computing real scores')
with open_file(out_fname + '.wig.gz', 'w') as out_signal_f:

    title = os.path.split(out_fname)[1].replace('.wig', '').replace('.gz', '')
    out_signal_f.write('track type=wiggle_0 name="%s" description="%s"\n' % (title, title))

    closest_transition_cache = {}
    scores_cache = {}
    pure_state_transitions = generate_pure_state_transitions(states, metric_A, metric_B)
    for chrom in chromosomes:
        echo('Chromosome:', chrom)

        signal, state_transitions = get_overall_and_per_state_diff_score( chrom,
                                                       group_A_segmentations,
                                                       metric_A,
                                                       group_B_segmentations,
                                                       metric_B,
                                                       BIN_SIZE,
                                                       states,
                                                       scores_cache,
                                                       pure_state_transitions,
                                                       closest_transition_cache)

        # signal, state_transitions = get_diff_score(chrom,
        #                                            group_A_segmentations,
        #                                            metric_A,
        #                                            group_B_segmentations,
        #                                            metric_B,
        #                                            BIN_SIZE,
        #                                            pure_state_transitions,
        #                                            closest_transition_cache,
        #                                            scores_cache)
        if args.smooth:
            signal = smooth(signal)

        real_scores[chrom] = signal
        real_state_transitions[chrom] = state_transitions

        store_wig(chrom,
                  signal,
                  out_signal_f=out_signal_f,
                  span=BIN_SIZE)


# determine FDR threshold
all_scores = sorted(s for c in real_scores for s in real_scores[c])

score_cutoff = None
best_score_cutoff = 0
best_fdr = 1

score_idx = 0
bgr_idx = 0

score_to_fdr = {}

while True:

    # if we reached the end of all_scores, then nothing is significant
    if score_idx == len(all_scores):
        break

    # find how many background scores are greater than the current real score
    while bgr_idx < len(background_scores) and background_scores[bgr_idx] < all_scores[score_idx]:
        bgr_idx += 1

    false_positives = (len(background_scores) - bgr_idx) / float(len(background_scores))
    all_positives = (len(all_scores) - score_idx + 1) / float(len(all_scores))
    current_fdr = false_positives / all_positives

    if current_fdr < best_fdr:
        best_fdr = current_fdr
        best_score_cutoff = all_scores[score_idx]
    else:
        current_fdr = best_fdr

    score_to_fdr[all_scores[score_idx]] = current_fdr

    if current_fdr <= fdr_threshold and score_cutoff is None:
        score_cutoff = all_scores[score_idx]

    score_idx += 1

    # skip identical scores
    while score_idx < len(all_scores) and all_scores[score_idx] == all_scores[score_idx - 1]:
        score_idx += 1


#
#
#
#
# score_cutoff = None
# score_idx = 0
# bgr_idx = 0
#
# while True:
#
#     # skip identical score values
#     if 0 < score_idx < len(all_scores) and all_scores[score_idx] == all_scores[score_idx - 1]:
#         score_idx += 1
#
#     # if we reached the end of all_scores, then nothing is significant
#     if score_idx == len(all_scores):
#         score_cutoff = None
#         break
#
#     # find how many background scores are greater than the current real score
#     while bgr_idx < len(background_scores) and background_scores[bgr_idx] < all_scores[score_idx]:
#         bgr_idx += 1
#
#     if bgr_idx < len(background_scores):
#         false_positives = (len(background_scores) - bgr_idx) / float(len(background_scores))
#         all_positives = (len(all_scores) - score_idx) / float(len(all_scores))
#
#         if false_positives / all_positives <= fdr_threshold:
#             score_cutoff = all_scores[score_idx]
#             break
#     else:
#         score_cutoff = all_scores[score_idx]
#         break
#
#     score_idx += 1

echo('Score cutoff at FDR', fdr_threshold, ':', score_cutoff)

if score_cutoff is None:
    echo('No significant changes were found at FDR of', fdr_threshold, '. Best FDR=', best_fdr, 'cutoff=', best_score_cutoff)
    score_cutoff = best_score_cutoff
    fdr_threshold = best_fdr

else:
    with open_file(out_fname + '.differential_peaks_FDR_' + str(fdr_threshold) + '.bed.gz', 'w') as out_diff_peaks_f:
        for chrom in sorted(real_scores):

            peak_start = None
            peak_no = 0
            peak_score = 0

            for bin_no, bin_signal in enumerate(real_scores[chrom]):
                bin_start = bin_no * BIN_SIZE
                if bin_signal >= score_cutoff:
                    if peak_start is None:
                        peak_start = bin_start
                        peak_score = 0
                        peak_no += 1

                    peak_score += bin_signal

                if (bin_signal < score_cutoff or bin_no == len(real_scores[chrom]) - 1) and peak_start is not None:

                    if bin_no == len(real_scores[chrom]) - 1:
                        bin_start += BIN_SIZE

                    out_diff_peaks_f.write('\t'.join(map(str, [chrom,
                                                               peak_start,
                                                               bin_start,
                                                               chrom + '_' + str(peak_no),
                                                               '%.2lf' % (BIN_SIZE * peak_score / (bin_start - peak_start)),
                                                               '.'])) + '\n')
                    peak_start = None

with open_file(out_fname + '.bed.gz', 'w') as out_states_f:
    for chrom in chromosomes:
        store_bed(chrom,
                  real_state_transitions[chrom],
                  real_scores[chrom],
                  score_to_fdr,
                  out_states_f)

echo('Storing scores in:', out_fname + '.scores.pickle')
with open(out_fname + '.scores.pickle', 'w') as scores_f:
    pickle.dump({'bgr': background_scores,
                 'fgr': all_scores}, scores_f, protocol=pickle.HIGHEST_PROTOCOL)
    # pickle.dump(all_scores, scores_f, protocol=pickle.HIGHEST_PROTOCOL)
#
# from matplotlib import pyplot as plt
# hist_info = plt.hist([background_scores, all_scores], bins=50, color=['blue', 'red'], label=['bgr', 'real'],
#                      normed=1)
# # max_y = max(v for a in [hist_info[0][0][2:], hist_info[0][1][2:]] for v in a)
#
# # xmin, xmax = plt.xlim()
# # plt.xlim(1, xmax)
# # plt.ylim(0, max_y)
# plt.legend()
# plt.savefig(out_fname + '_score_histogram.png', dpi=300)
#
echo('Done')
