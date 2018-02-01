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
                                       posteriors_dir):
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

        # compute the score between A and B
        score.append(symmetric_KL_divergence(a_probs, b_probs))
        state_transitions.append(state_transition(a_probs, b_probs))
        # print 'A:', a_bins
        # print 'B:', b_bins
        # print

    for f in a_posteriors_files + b_posteriors_files:
        f.close()

    return score, state_transitions


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
                              to_smooth):

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

    background_scores = []
    seen = set([(A_segs_keys, B_segs_keys),
                (B_segs_keys, A_segs_keys)])
    
    for shuffled_A_segs, shuffled_B_segs in shuffled_samples():

        if (shuffled_A_segs, shuffled_B_segs) in seen or (shuffled_B_segs, shuffled_A_segs) in seen:
            continue

        seen.add((shuffled_A_segs, shuffled_B_segs))
        seen.add((shuffled_B_segs, shuffled_A_segs))

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

        signal, _ = get_diff_score_based_on_posteriors(longest_chromosome,
                                                       reduced_A_segs,
                                                       metric_A,
                                                       reduced_B_segs,
                                                       metric_B,
                                                       BIN_SIZE,
                                                       posteriors_dir)
        if to_smooth:
            signal = smooth(signal)

        background_scores.extend(signal)

    return sorted(background_scores)

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
print_consistent_state_transition_scores(states, metric_A, metric_B)

fdr_threshold = args.fdr_threshold


if fdr_threshold < 1:
    echo('Learning significance threshold for at p-value:', fdr_threshold)
    background_scores = compute_background_scores(group_A_segmentations,
                                                  group_B_segmentations,
                                                  chromosomes,
                                                  states,
                                                  posteriors_dir,
                                                  args.smooth)

# longest_chromosome = max(chromosomes, key=lambda c: group_A_segmentations.values()[0][c][-1][1])
# print 'significance threshold at', fdr_threshold, ':', score_cutoff
# exit(1)

out_fname = args.o


# threshold = determine_threshold(group_A_segmentations, metric_A,
#                                 group_B_segmentations, metric_B,
#                                 fdr_threshold=.01)
real_scores = {}
echo('Storing output:', out_fname)
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
                                                                       posteriors_dir)
        if args.smooth:
            signal = smooth(signal)

        real_scores[chrom] = signal

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

# determine FDR threshold
all_scores = sorted(s for c in real_scores for s in real_scores[c])
score_cutoff = None
score_idx = 0
bgr_idx = 0

while True:

    # skip identical score values
    if 0 < score_idx < len(all_scores) and all_scores[score_idx] == all_scores[score_idx - 1]:
        score_idx += 1

    # if we reached the end of all_scores, then nothing is significant
    if score_idx == len(all_scores):
        score_cutoff = None
        break

    # find how many background scores are greater than the current real score
    while bgr_idx < len(background_scores) and background_scores[bgr_idx] < all_scores[score_idx]:
        bgr_idx += 1

    if bgr_idx < len(background_scores):
        false_positives = (len(background_scores) - bgr_idx) / float(len(background_scores))
        all_positives = (len(all_scores) - score_idx) / float(len(all_scores))

        if false_positives / all_positives <= fdr_threshold:
            score_cutoff = all_scores[score_idx]
            break
    else:
        score_cutoff = all_scores[score_idx]
        break

    score_idx += 1

echo('Score cutoff at FDR', fdr_threshold, ':', score_cutoff)
significant_peaks = []

if score_cutoff is None:
    echo('No significant changes were found!')
else:
    with open_file(out_fname + '.differential_peaks.bed.gz', 'w') as out_diff_peaks_f:
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

                    peak_score = (BIN_SIZE * peak_score / (bin_start - peak_start))
                    peak_id = chrom + '_' + str(peak_no)
                    significant_peaks.append([chrom, peak_start, bin_start, peak_id, peak_score])

                    out_diff_peaks_f.write('\t'.join(map(str, [chrom,
                                                               peak_start,
                                                               bin_start,
                                                               peak_id,
                                                               '%.2lf' % peak_score,
                                                               '.'])) + '\n')
                    peak_start = None


echo('Cluster significant peaks')
K_clusters = 10
MAX_K_MEANS_ITERATIONS = 200
k_means_iteration = 0
total_distance = 0
groups = ['A', 'B']
priors = None

get_new_cluster = lambda is_random: dict((group,
                                            dict((state, random.random() if is_random else 0) for state in states))
                                           for group in groups)


def nearest_cluster(clusters, state_frequencies):
    min_cluster_idx = -1
    min_cluster_distance = float('inf')

    for cluster_idx, cluster in enumerate(clusters):
        _flat_dict = lambda d: dict((g + ' ' + k, d[g][k]) for g in d for k in d[g])
        # try:
        # dist = hellinger_distance(_flat_dict(state_frequencies), _flat_dict(cluster))
        dist = KL(_flat_dict(state_frequencies), _flat_dict(cluster))
        # except:
        #     print _flat_dict(state_frequencies)
        #     print _flat_dict(cluster)
        #     raise

        if dist <= min_cluster_distance:
            min_cluster_distance = dist
            min_cluster_idx = cluster_idx

    return min_cluster_idx, min_cluster_distance

cur_chrom = None
cur_chrom_position = 0
cur_chrom_bins_iterator = None
significant_bins = []



datasets_A = sorted(group_A_segmentations)
datasets_B = sorted(group_B_segmentations)


for peak_idx in xrange(len(significant_peaks)):
    state_frequencies = get_new_cluster(is_random=False)

    chrom, peak_start, peak_end, peak_id, peak_score = significant_peaks[peak_idx]

    if chrom != cur_chrom:
        cur_chrom = chrom
        cur_chrom_position = 0
        # cur_chrom_bins_iterator = izip(izip(*[iterstate(group_A_segmentations[d][chrom], BIN_SIZE) for d in group_A_segmentations]),
        #                                izip(*[iterstate(group_B_segmentations[d][chrom], BIN_SIZE) for d in group_B_segmentations]))
        #
        a_posteriors_files = [iter_posteriors(get_chrom_posteriors_fname(posteriors_dir, d, chrom)) for d in datasets_A]
        b_posteriors_files = [iter_posteriors(get_chrom_posteriors_fname(posteriors_dir, d, chrom)) for d in datasets_B]

        cur_chrom_bins_iterator = izip(izip(*a_posteriors_files),
                                       izip(*b_posteriors_files))

    while cur_chrom_position < peak_start:

        a_posteriors, b_posteriors = cur_chrom_bins_iterator.next()
        cur_chrom_position += BIN_SIZE

    while cur_chrom_position < peak_end:
        a_posteriors, b_posteriors = cur_chrom_bins_iterator.next()
        a_probs = get_prob_vector_based_on_posteriors(a_posteriors, datasets_A, metric_A)
        b_probs = get_prob_vector_based_on_posteriors(b_posteriors, datasets_B, metric_B)
        significant_bins.append([chrom,
                                 cur_chrom_position,
                                 cur_chrom_position + BIN_SIZE,
                                 peak_id + '_' + str((cur_chrom_position - peak_start) / BIN_SIZE + 1),
                                 symmetric_KL_divergence(a_probs, b_probs),
                                 {'A': a_probs, 'B': b_probs}])

        cur_chrom_position += BIN_SIZE


clusters = map(copy.deepcopy,
               random.sample([state_frequencies for _, _, _, _, _, state_frequencies in significant_bins], K_clusters))

prev_total_distance = 0
while k_means_iteration < MAX_K_MEANS_ITERATIONS:
    k_means_iteration += 1
    print 'Iteration:', k_means_iteration

    total_distance = 0
    delta = 0
    cluster_means = [get_new_cluster(is_random=False) for _ in xrange(K_clusters)]

    cluster_elements = [0] * K_clusters
    cluster_scores = [0] * K_clusters

    for chrom, peak_start, peak_end, peak_id, peak_score, state_frequencies in significant_bins:

        min_cluster_idx, min_cluster_distance = nearest_cluster(clusters, state_frequencies)
        total_distance += min_cluster_distance
        cluster_elements[min_cluster_idx] += 1
        cluster_scores[min_cluster_idx] += peak_score

        for group in ['A', 'B']:
            for state in states:
                cluster_means[min_cluster_idx][group][state] += state_frequencies[group][state]

    for cluster_idx in xrange(K_clusters):
        if cluster_elements[cluster_idx] == 0:
            print cluster_idx, 'has no elements'
            continue
        for group in groups:
            for state in states:
                new_value = cluster_means[cluster_idx][group][state] / cluster_elements[cluster_idx]
                delta += abs(clusters[cluster_idx][group][state] - new_value)
                clusters[cluster_idx][group][state] = new_value

    priors = [float(n_elements) / len(significant_bins) for n_elements in cluster_elements]
    cluster_scores = [float(cluster_scores[cluster_idx]) / cluster_elements[cluster_idx]
                      for cluster_idx in xrange(K_clusters)]

    print 'Delta:', delta
    print 'Total distance:', total_distance, 'delta:', prev_total_distance - total_distance
    prev_total_distance = total_distance
    for cluster_idx in xrange(len(clusters)):
        print '\t'.join(map(str, [cluster_idx, priors[cluster_idx], cluster_scores[cluster_idx]]))
        print '\t'.join(['State'] + states)
        print '\t'.join(['A'] + ['%.2lf' % v for s, v in sorted(clusters[cluster_idx]['A'].items(), key=lambda (s, v): state_key(s))])
        print '\t'.join(['B'] + ['%.2lf' % v for s, v in sorted(clusters[cluster_idx]['B'].items(), key=lambda (s, v): state_key(s))])
    # pprint.pformat(clusters)

    echo('Priors:', priors)

    if delta < 0.1 or k_means_iteration == MAX_K_MEANS_ITERATIONS:
        cluster_files = [open_file(out_fname + '_cluster_' + str(cluster_idx) + '.bed.gz', 'w') for cluster_idx in xrange(K_clusters)]

        for chrom, peak_start, peak_end, peak_id, peak_score, state_frequencies in significant_bins:
            min_cluster_idx, min_cluster_distance = nearest_cluster(clusters, state_frequencies)
            cluster_files[min_cluster_idx].write('\t'.join(map(str, [chrom, peak_start, peak_end, peak_id, peak_score, '.'])) + '\n')




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
