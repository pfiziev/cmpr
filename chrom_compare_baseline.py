import argparse
from itertools import izip, combinations
from multiprocessing import Pool
import os
import pprint
import random
import sys
import operator
import pickle
import itertools
import gc
import traceback

from utils import *
from metric import *


def _map(func, array):
    for arg in array:
        yield func(arg)


def compress_segmentations_into_non_repeating_windows(chrom_segmentations_A, chrom_segmentations_B, max_min_window):

    chrom_length = len(chrom_segmentations_A[0])

    non_repating_windows = []
    prev_a_windows = None
    prev_b_windows = None
    prev_a_bins = None
    prev_b_bins = None

    unique_windows = set()
    unique_windows_indices = {}

    for bin_idx in xrange(chrom_length):

        a_bins = tuple(seg[bin_idx] for seg in chrom_segmentations_A)
        b_bins = tuple(seg[bin_idx] for seg in chrom_segmentations_B)

        window_start = max(0, bin_idx - max_min_window)
        window_end = min(chrom_length, bin_idx + max_min_window + 1)

        a_window = tuple(s for c in sorted(set(tuple(seg[i] for seg in chrom_segmentations_A) for i in xrange(window_start, window_end))) for s in c)
        b_window = tuple(s for c in sorted(set(tuple(seg[i] for seg in chrom_segmentations_B) for i in xrange(window_start, window_end))) for s in c)
        # b_window = tuple(seg[i] for i in xrange(window_start, window_end) for seg in chrom_segmentations_B)

        # a_window = tuple(seg[i] for i in xrange(window_start, window_end) for seg in chrom_segmentations_A)
        # b_window = tuple(seg[i] for i in xrange(window_start, window_end) for seg in chrom_segmentations_B)

        if a_window != prev_a_windows or b_window != prev_b_windows or a_bins != prev_a_bins or b_bins != prev_b_bins:

            key = (a_bins, b_bins, a_window, b_window)

            if key not in unique_windows:
                unique_windows_indices[key] = len(unique_windows)
                unique_windows.add(key)

            non_repating_windows.append([unique_windows_indices[key], 0])

            prev_a_bins = a_bins
            prev_b_bins = b_bins
            prev_a_windows = a_window
            prev_b_windows = b_window


        non_repating_windows[-1][-1] += 1

    return non_repating_windows, sorted(unique_windows, key=lambda k: unique_windows_indices[k])


def expand_non_overlapping_windows(chrom, scores, state_transitions, chrom_segmentations_windows, background_chunk=False):

    total_chrom_length = sum(c for _, c in chrom_segmentations_windows)
    score_types = scores[0].keys()

    expanded_scores = dict((score_type, [0.0] * total_chrom_length) for score_type in score_types)
    expanded_state_transitions = None
    state_transition = None

    if not background_chunk:
        expanded_state_transitions = [None] * total_chrom_length
    bin_idx = 0

    for win_idx in xrange(len(chrom_segmentations_windows)):

        if not background_chunk:
            state_transition = state_transitions[win_idx]

        _, count = chrom_segmentations_windows[win_idx]

        for _ in xrange(count):
            if not background_chunk:
                expanded_state_transitions[bin_idx] = state_transition

            for score_type in score_types:
                expanded_scores[score_type][bin_idx] = scores[win_idx][score_type]

            bin_idx += 1

    return expanded_scores if background_chunk else (chrom, expanded_scores, expanded_state_transitions)

MAX_CACHE_SIZE = 10000


def get_prob_vector_marks(bins, datasets, metric, states):

    probs = dict((s, 0) for s in states)

    for s in states:
        for dataset, d_state in izip(datasets, bins):
            probs[s] += metric[dataset][d_state][s]

    total = len(datasets)

    for s in states:
        probs[s] /= total

    return probs


def read_ChromHMM_posteriors(chrom, posteriors_dir, datasets):
    celltypes = sorted([os.path.split(d)[1].split('_segments')[0] for d in datasets])

    signal = {}

    for ct in celltypes:
        echo('Reading posterior for:', ct, chrom)
        fname = os.path.join(posteriors_dir, ct + '_' + chrom + '_posterior.txt.gz')

        with open_file(fname) as in_f:
            n_bins = len(list(in_f)) - 2

        with open_file(fname) as in_f:
            ct, chrom = in_f.readline().strip().split()

            if ct not in signal:
                signal[ct] = {}

            states = in_f.readline().strip().split()

            for s in states:
                signal[ct][s] = [0] * n_bins

            for i, l in enumerate(in_f):
                for s, v in izip(states, map(float, l.strip().split())):
                    signal[ct][s][i] = round(v, 3)
                if i < 2:
                    print dict((s, signal[ct][s][i]) for s in states)

    return signal


def compute_posteriors_diff(chrom_posteriors_A, chrom_posteriors_B, compute_per_state_scores, use_symmetric_KL=True):
    n_bins = len(chrom_posteriors_A.values()[0].values()[0])
    states = sorted(chrom_posteriors_A.values()[0])
    signal = dict((score_type, [0] * n_bins) for score_type in [OVERALL_SCORE] + (states if compute_per_state_scores else []))

    for bin_i in xrange(n_bins):
        p_A = dict((s, mean([chrom_posteriors_A[ct][s][bin_i] for ct in chrom_posteriors_A])) for s in states)
        p_B = dict((s, mean([chrom_posteriors_B[ct][s][bin_i] for ct in chrom_posteriors_B])) for s in states)

        def pseudo(p):
            PSEUDO = 0.001
            total = sum(p.values()) + PSEUDO * len(p)
            return dict((s, (p[s] + PSEUDO) / total) for s in p)

        # if bin_i < 2:
        #     print pseudo(p_A)
        #     print pseudo(p_B)
        #     print

        # signal[OVERALL_SCORE][bin_i] = hellinger_distance(p_A, p_B)
        if use_symmetric_KL:
            signal[OVERALL_SCORE][bin_i] = symmetric_KL_divergence(pseudo(p_A), pseudo(p_B))
        else:
            signal[OVERALL_SCORE][bin_i] = hellinger_distance(p_A, p_B)

        if compute_per_state_scores:
            for s in states:
                signal[s][bin_i] = p_A[s] - p_B[s]

    return signal


def get_overall_and_per_state_diff_score_posteriors(chrom,
                                                    group_A_segmentations,
                                                    group_B_segmentations,
                                                    posteriors_dir,
                                                    BIN_SIZE,
                                                    compute_per_state_scores,
                                                    use_symmetric_KL):
    datasets_A = sorted(group_A_segmentations)
    datasets_B = sorted(group_B_segmentations)

    chrom_posteriors_A = read_ChromHMM_posteriors(chrom, posteriors_dir, datasets_A)
    chrom_posteriors_B = read_ChromHMM_posteriors(chrom, posteriors_dir, datasets_B)
    signal = compute_posteriors_diff(chrom_posteriors_A, chrom_posteriors_B, compute_per_state_scores, use_symmetric_KL)
    state_transitions = []
    for (a_bins, b_bins) in izip(izip(*[iterstate(group_A_segmentations[d][chrom], BIN_SIZE) for d in datasets_A]),
                                 izip(*[iterstate(group_B_segmentations[d][chrom], BIN_SIZE) for d in datasets_B])):
        state_transitions.append(('X', 'X', a_bins, b_bins))

    return chrom, signal, state_transitions


def get_overall_and_per_state_diff_score_by_state_counts(chrom,
                                                         group_A_segmentations,
                                                         group_B_segmentations,
                                                         BIN_SIZE,
                                                         states,
                                                         compute_per_state_scores):

    echo('Processing chromosome:', chrom)

    datasets_A = sorted(group_A_segmentations)
    datasets_B = sorted(group_B_segmentations)

    n_bins = (group_A_segmentations.values()[0][chrom][-1][1]) / BIN_SIZE

    signal = dict((score_type, [0] * n_bins) for score_type in [OVERALL_SCORE] + (states if compute_per_state_scores else []))

    state_transitions = []
    for bin_idx, (a_bins, b_bins) in enumerate(izip(izip(*[iterstate(group_A_segmentations[d][chrom], BIN_SIZE) for d in datasets_A]),
                                                    izip(*[iterstate(group_B_segmentations[d][chrom], BIN_SIZE) for d in datasets_B]))):
        a_set = set(a_bins)
        b_set = set(b_bins)

        union = a_set | b_set

        signal[OVERALL_SCORE][bin_idx] = float(len([s for s in a_bins if s not in b_set])) / len(a_bins) + float(len([s for s in b_bins if s not in a_set])) / len(b_bins)

        if compute_per_state_scores:
            for state in union:
                score = float(len([s for s in a_bins if s == state])) / len(a_bins) - float(len([s for s in b_bins if s == state])) / len(b_bins)
                signal[state][bin_idx] = score

        state_transitions.append(('X', 'X', a_bins, b_bins))

    return chrom, signal, state_transitions


def get_overall_and_per_state_diff_score_by_fishers_exact_test(chrom,
                                                               group_A_segmentations,
                                                               group_B_segmentations,
                                                               BIN_SIZE,
                                                               states,
                                                               compute_per_state_scores,
                                                               keep_max_scores=False):

    echo('Processing chromosome:', chrom)
    cache = dict((s, {}) for s in states)

    from scipy.stats import fisher_exact
    datasets_A = sorted(group_A_segmentations)
    datasets_B = sorted(group_B_segmentations)

    n_bins = (group_A_segmentations.values()[0][chrom][-1][1]) / BIN_SIZE

    signal = dict((score_type, [0] * n_bins) for score_type in [OVERALL_SCORE] + (states if compute_per_state_scores else []))

    state_transitions = []
    for bin_idx, (a_bins, b_bins) in enumerate(izip(izip(*[iterstate(group_A_segmentations[d][chrom], BIN_SIZE) for d in datasets_A]),
                                                    izip(*[iterstate(group_B_segmentations[d][chrom], BIN_SIZE) for d in datasets_B]))):
        a_set = set(a_bins)
        b_set = set(b_bins)
        union = a_set | b_set

        signal[OVERALL_SCORE][bin_idx] = float(len([s for s in a_bins if s not in b_set])) / len(a_bins) + float(len([s for s in b_bins if s not in a_set])) / len(b_bins)

        if compute_per_state_scores:
            for state in union:
                a = len([s for s in a_bins if s == state])
                b = len(a_bins) - a
                c = len([s for s in b_bins if s == state])
                d = len(b_bins) - c

                key = (a, c)

                if key not in cache[state]:
                    odds_ratio, pval = fisher_exact([[a, b], [c, d]])
                    score = -math.log(pval, 2) * sign(float(a) / (a + b) - float(c) / (c + d))

                    if score == 0:
                        score = 0.00001

                    cache[state][key] = score

                signal[state][bin_idx] = cache[state][key]

            if keep_max_scores:
                min_score_state = min(union, key=lambda s: signal[s][bin_idx])
                max_score_state = max(union, key=lambda s: signal[s][bin_idx])

                for state in union:
                    if state not in [max_score_state, min_score_state]:
                        signal[state][bin_idx] = 0

        state_transitions.append(('X', 'X', a_bins, b_bins))

    return chrom, signal, state_transitions



def get_overall_and_per_state_diff_score_signal( chrom,
                                                 group_A_segmentations,
                                                 group_B_segmentations,
                                                 signal_dir,
                                                 BIN_SIZE,
                                                 LOG2RPKM,
                                                 total_reads_per_mark,
                                                 compute_per_state_scores,
                                                 max_min_window):

    datasets_A = sorted(group_A_segmentations)
    datasets_B = sorted(group_B_segmentations)

    chrom_signal_A = read_ChromHMM_signal(chrom, signal_dir, datasets_A, total_reads_per_mark, BIN_SIZE, LOG2RPKM)
    chrom_signal_B = read_ChromHMM_signal(chrom, signal_dir, datasets_B, total_reads_per_mark, BIN_SIZE, LOG2RPKM)
    signal = compute_signal_diff(chrom_signal_A, chrom_signal_B, compute_per_state_scores=compute_per_state_scores)
    state_transitions = []
    for (a_bins, b_bins) in izip(izip(*[iterstate(group_A_segmentations[d][chrom], BIN_SIZE) for d in datasets_A]),
                                 izip(*[iterstate(group_B_segmentations[d][chrom], BIN_SIZE) for d in datasets_B])):
        state_transitions.append(('X', 'X', a_bins, b_bins))

    return chrom, signal, state_transitions

def get_overall_and_per_state_diff_score_emissions(chrom,
                                         group_A_segmentations,
                                         group_B_segmentations,
                                         emissions,
                                         BIN_SIZE,
                                         states,
                                         compute_per_state_scores,
                                         max_min_window,
                                         background_chunk=False):

    if not background_chunk:
        echo('Chromosome:', chrom)

    cache_miss = 0

    datasets_A = sorted(group_A_segmentations)
    datasets_B = sorted(group_B_segmentations)

    if background_chunk:
        chrom_segmentations_A = [group_A_segmentations[d] for d in datasets_A]
        chrom_segmentations_B = [group_B_segmentations[d] for d in datasets_B]
    else:
        chrom_segmentations_A = [chrom_segmentation_to_list(group_A_segmentations[d][chrom], BIN_SIZE) for d in datasets_A]
        chrom_segmentations_B = [chrom_segmentation_to_list(group_B_segmentations[d][chrom], BIN_SIZE) for d in datasets_B]

    if not background_chunk:
        echo('Chrom length:', len(chrom_segmentations_A[0]))

    chrom_segmentations_windows, unique_windows = compress_segmentations_into_non_repeating_windows(
                                                                                    chrom_segmentations_A,
                                                                                    chrom_segmentations_B,
                                                                                    max_min_window)

    del chrom_segmentations_A
    del chrom_segmentations_B

    if not background_chunk:
        gc.collect()
        echo('Done compressing into windows:', len(chrom_segmentations_windows), 'unique windows:', len(unique_windows))

    chrom_length = len(chrom_segmentations_windows)
    scores = [None] * chrom_length # dict((k, [0.0] * chrom_length) for k in [OVERALL_SCORE] + (states if compute_per_state_scores else []))

    scores_cache = [None] * len(unique_windows)
    closest_transition_cache = [None] * len(unique_windows)

    state_transitions = None
    if not background_chunk:
        state_transitions = [None] * chrom_length

    # return scores, state_transitions
    marks = sorted(emissions.values()[0])
    metric_A = dict((d, emissions) for d in datasets_A)
    metric_B = dict((d, emissions) for d in datasets_B)

    for bin_idx, (window_idx, _) in enumerate(chrom_segmentations_windows):

        a_bins, b_bins, a_window, b_window = unique_windows[window_idx]

        if scores_cache[window_idx] is None:
            cache_miss += 1

            # get the posterior probabilities for group A by summing up the
            # the probability vectors for each state

            a_probs = get_prob_vector_marks(a_bins, datasets_A, metric_A, marks)

            # same for group B
            b_probs = get_prob_vector_marks(b_bins, datasets_B, metric_B, marks)
            score = eucld([a_probs[m] for m in marks], [b_probs[m] for m in marks])

            scores_dict = dict((k, 0) for k in [OVERALL_SCORE] + (states if compute_per_state_scores else []))
            scores_dict[OVERALL_SCORE] = score

            # if compute_per_state_scores:
            #     for state in states:
            #         scores_dict[state] = true_a_probs[state] - true_b_probs[state]

            if not background_chunk:
                # flat_prob_dict = dict((group + ' ' + key, d[key])
                #                       for group, d in [('A', true_a_probs), ('B', true_b_probs)] for key in d)
                #
                # true_a_bins = set(a_bins if true_a_idx == -1 else a_window[true_a_idx * len(datasets_A):(true_a_idx + 1) * len(datasets_A)])
                # true_b_bins = set(b_bins if true_b_idx == -1 else b_window[true_b_idx * len(datasets_B):(true_b_idx + 1) * len(datasets_B)])
                #
                # state_pairs = [(s1, s2) for s1 in true_a_bins for s2 in true_b_bins]
                #
                # closest_state_A, closest_state_B = min(state_pairs,
                #                                        key=lambda (_s1, _s2):
                #                                        KL(flat_prob_dict,
                #                                           pure_state_transitions[_s1][_s2][STATE_DISTRIBUTION]))

                closest_transition_cache[window_idx] = ('X', 'X', a_bins, b_bins)

            scores_cache[window_idx] = scores_dict

        else:
            scores_dict = scores_cache[window_idx]

        # compute the score between A and B
        scores[bin_idx] = scores_dict

        if not background_chunk:
            state_transitions[bin_idx] = closest_transition_cache[window_idx]

    if not background_chunk:
        echo(chrom_length, len(scores_cache), 100 * (cache_miss / float(chrom_length)), cache_miss)

    return expand_non_overlapping_windows(chrom, scores, state_transitions, chrom_segmentations_windows, background_chunk)


def store_wig(chrom,
              signal,
              out_signal_f,
              span):

    out_signal_f.write('fixedStep\tchrom=%s\tstart=0\tstep=%d\tspan=%d\n' % (chrom, span, span))

    for bin_signal in signal:

        out_signal_f.write(str(bin_signal) + '\n')


def store_dict_wig(chrom,
                   signal,
                   out_signal_files,
                   span):

    for signal_type in signal:
        # echo(signal_type)
        store_wig(chrom, signal[signal_type], out_signal_files[signal_type], span)


def store_dict_bed(chrom,
                   state_transitions,
                   signal,
                   background_model,
                   out_summary_files,
                   span):

    for signal_type in sorted(signal):
        # echo(signal_type)
        store_bed(chrom,
                  state_transitions,
                  signal[signal_type],
                  out_summary_files[signal_type],
                  background_model,
                  signal_type,
                  span)

import bisect

get_FDR_cache = {}
def get_FDR(signal, background_model, signal_type):
    key = (signal, signal_type)
    if key not in get_FDR_cache:
        bm = background_model[signal_type]
        idx = bisect.bisect_left(bm, (signal, 0))

        if idx == len(bm):
            get_FDR_cache[key] = bm[-1][1]
        else:
            bm_signal, bm_fdr = bm[idx]

            if signal == bm_signal:
                get_FDR_cache[key] = bm_fdr
            else:
                # interpolate the FDR
                prev_signal, prev_fdr = bm[idx - 1]
                get_FDR_cache[key] = ((bm_fdr - prev_fdr) / (bm_signal - prev_signal)) * signal + (prev_fdr * bm_signal - bm_fdr * prev_signal) / (bm_signal - prev_signal)

    return get_FDR_cache[key]



def store_bed(chrom,
              state_transitions,
              signal,
              out_summary_file,
              background_model,
              signal_type,
              span):

    start = 0
    for bin_no, (closest_state_A,
                 closest_state_B,
                 a_states,
                 b_states) in enumerate(state_transitions):

        bin_signal = signal[bin_no]
        out_summary_file.write('\t'.join(map(str,
                                         [chrom,
                                          start,
                                          start + span,
                                          '.',
                                          bin_signal,
                                          '.',
                                          closest_state_A + '-' + closest_state_B,
                                          ','.join(a_states) + '/' + ','.join(b_states),
                                          get_FDR(abs(bin_signal), background_model, signal_type)])) + '\n')

        start += span


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
                                                         random_chunks,
                                                         states,
                                                         emissions=None,
                                                         signal_dir=None,
                                                         total_reads_per_mark=None,
                                                         n_perm=100,
                                                         to_smooth=False,
                                                         max_min_window=0,
                                                         compute_per_state_scores=False,
                                                         n_threads=1,
                                                         BIN_SIZE=None,
                                                         LOG2RPKM=None
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

    if n_A_segs == n_B_segs:
        n_comb = n_combinations(len(all_seg_keys), n_A_segs) / 2 - 1
    else:
        n_comb = n_combinations(len(all_seg_keys), n_A_segs) - 1

    def shuffled_samples(n_perm):

        seen = set()

        if n_comb > n_perm:

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

    background_scores = dict((k, {}) for k in [OVERALL_SCORE] + (states if compute_per_state_scores else []))

    if n_threads == 1:
        _map = itertools.imap
    else:
        pool = Pool(processes=n_threads)
        _map = pool.imap_unordered

    for combo_background_scores in _map(worker,
                                        [(process_shuffled_combo,
                                          combo_no,
                                          shuffled_A_datasets,
                                          dict((dataset, all_segs[dataset]) for dataset in shuffled_A_datasets),

                                          shuffled_B_datasets,
                                          dict((dataset, all_segs[dataset]) for dataset in shuffled_B_datasets),

                                          emissions,
                                          signal_dir,
                                          total_reads_per_mark,

                                          random_chunks,
                                          BIN_SIZE,
                                          LOG2RPKM,
                                          states,
                                          compute_per_state_scores,
                                          max_min_window,
                                          to_smooth)
                                         for combo_no, (shuffled_A_datasets, shuffled_B_datasets) in enumerate(shuffled_samples(n_perm))]):

        for score_type in combo_background_scores:
            for s in combo_background_scores[score_type]:
                if s not in background_scores[score_type]:
                    background_scores[score_type][s] = 0
                background_scores[score_type][s] += combo_background_scores[score_type][s]

    if n_threads > 1:
        pool.close()

    for score_type in background_scores:
        total_regions = sum(background_scores[score_type].itervalues())
        total_so_far = 0

        for s in sorted(background_scores[score_type], reverse=True):
            total_so_far += background_scores[score_type][s]
            background_scores[score_type][s] = total_so_far / float(total_regions)

    return dict((score_type, sorted(background_scores[score_type].items())) for score_type in background_scores)


def worker(args):
    func = None

    try:
        func = args[0]
        return func(*args[1:])

    except Exception, e:
        print 'Caught exception in output worker thread (combo_on: %s):' % str(args[1])
        print func

        echo(e)
        if hasattr(open_log, 'logfile'):
            traceback.print_exc(file=open_log.logfile)
        traceback.print_exc()

        print
        raise e


def read_ChromHMM_signal(chrom, signal_dir, datasets, total_reads_per_mark, BIN_SIZE, LOG2RPKM):
    echo('LOG2RKM:', LOG2RPKM)
    celltypes = sorted([get_celltype(d) for d in datasets])
    # celltypes = [re.sub(r'_(\d+)_segments.bed(.gz)?', '', os.path.split(d)[1]) for d in datasets]
    signal = {}

    for ct in celltypes:
        echo('Reading signal for:', ct, chrom)
        fname = os.path.join(signal_dir, ct + '_' + chrom + '_signal.txt.gz')

        with open_file(fname) as in_f:
            n_bins = len(list(in_f)) - 2

        with open_file(fname) as in_f:
            ct, chrom = in_f.readline().strip().split()

            if ct not in signal:
                signal[ct] = {}

            # if chrom not in signal[ct]:
            #     signal[ct][chrom] = {}

            marks = in_f.readline().strip().split()

            for m in marks:
                # signal[ct][chrom][m] = [0] * n_bins
                signal[ct][m] = [0] * n_bins

            for i, l in enumerate(in_f):
                for m, v in izip(marks, map(int, l.strip().split())):
                    # signal[ct][chrom][m][i] = v

                    signal[ct][m][i] = float(v * (10 ** 9 / BIN_SIZE)) / total_reads_per_mark[ct][m]

                    if LOG2RPKM:
                        signal[ct][m][i] = math.log(signal[ct][m][i] + 1, 2)

    return signal


def compute_signal_diff(chrom_signal_A, chrom_signal_B, chunk_start=None, chunk_end=None, compute_per_state_scores=False):

    start = 0
    end = len(chrom_signal_A.values()[0].values()[0])
    marks = sorted(chrom_signal_A.values()[0].keys())

    if chunk_start is not None:
        start = chunk_start
        end = chunk_end

    all_scores = dict((score_type, [0] * (end - start)) for score_type in [OVERALL_SCORE] + (marks if compute_per_state_scores else []))
    # scores = [0] * (end - start)

    for score_i, bin_i in enumerate(xrange(start, end)):
        all_scores[OVERALL_SCORE][score_i] = eucld([mean([chrom_signal_A[ct][m][bin_i] for ct in chrom_signal_A]) for m in marks],
                                                   [mean([chrom_signal_B[ct][m][bin_i] for ct in chrom_signal_B]) for m in marks])
        if compute_per_state_scores:
            for m in marks:
                all_scores[m][score_i] = mean([chrom_signal_A[ct][m][bin_i] for ct in chrom_signal_A]) - mean([chrom_signal_B[ct][m][bin_i] for ct in chrom_signal_B])

    return all_scores


def process_shuffled_combo( combo_no,
                            shuffled_A_datasets,
                            shuffled_A_segmentations,


                            shuffled_B_datasets,
                            shuffled_B_segmentations,

                            emissions,
                            signal_dir,
                            total_reads_per_mark,

                            random_chunks,
                            BIN_SIZE,
                            LOG2RPKM,
                            states,
                            compute_per_state_scores,
                            max_min_window,
                            to_smooth):

    echo('Combo no:', combo_no)

    print 'Shuffled A:', shuffled_A_datasets
    print 'Shuffled B:', shuffled_B_datasets
    print

    # shuffled_metric_A = learn_metric(shuffled_A_segmentations, states, BIN_SIZE)
    # print_metric(shuffled_metric_A)
    #
    # shuffled_metric_B = learn_metric(shuffled_B_segmentations, states, BIN_SIZE)

    # print_metric(shuffled_metric_B)
    # print '*' * 50
    # print_average_metric(states, shuffled_metric_A, shuffled_metric_B)
    # print '*' * 50

    # iterate over the bins in both segmentation groups
    # for chrom in random.sample(chromosomes, 2):

    background_scores = dict((k, {}) for k in [OVERALL_SCORE] + (states if compute_per_state_scores else []))

    for chrom in sorted(random_chunks):
        # echo(combo_no, chrom)
        if emissions:
            chrom_segmentations_A = [chrom_segmentation_to_list(shuffled_A_segmentations[d][chrom], BIN_SIZE) for d in shuffled_A_datasets]
            chrom_segmentations_B = [chrom_segmentation_to_list(shuffled_B_segmentations[d][chrom], BIN_SIZE) for d in shuffled_B_datasets]
        elif signal_dir:
            chrom_signal_A = read_ChromHMM_signal(chrom, signal_dir, shuffled_A_datasets, total_reads_per_mark, BIN_SIZE, LOG2RPKM)
            chrom_signal_B = read_ChromHMM_signal(chrom, signal_dir, shuffled_B_datasets, total_reads_per_mark, BIN_SIZE, LOG2RPKM)

        for chunk_start, chunk_end in random_chunks[chrom]:
            # echo(combo_no, chunk_start, chunk_end)
            chunk_scores = None
            if emissions:
                chunk_segmentations_A = dict((d, seg[chunk_start:chunk_end]) for d, seg in izip(shuffled_A_datasets, chrom_segmentations_A))
                chunk_segmentations_B = dict((d, seg[chunk_start:chunk_end]) for d, seg in izip(shuffled_B_datasets, chrom_segmentations_B))

                chunk_scores = get_overall_and_per_state_diff_score_emissions(chrom,
                                                                    chunk_segmentations_A,
                                                                    chunk_segmentations_B,
                                                                    emissions,
                                                                    BIN_SIZE,
                                                                    states,
                                                                    compute_per_state_scores,
                                                                    max_min_window,
                                                                    background_chunk=True)
            elif signal_dir:
                chunk_scores = compute_signal_diff(chrom_signal_A, chrom_signal_B, chunk_start, chunk_end, compute_per_state_scores=compute_per_state_scores)

            if to_smooth:
                chunk_scores = dict((score_type, smooth(chunk_scores[score_type])) for score_type in chunk_scores)

            for score_type in background_scores:
                for s in chunk_scores[score_type]:
                    abs_score = abs(s)
                    if abs_score not in background_scores[score_type]:
                        background_scores[score_type][abs_score] = 0

                    background_scores[score_type][abs_score] += 1
    gc.collect()
    return background_scores

def compute_foreground_scores(group_A_segmentations,
                              group_B_segmentations,
                              random_chunks,
                              states,

                              emissions,
                              signal_dir,
                              total_reads_per_mark,
                              BIN_SIZE=None,
                              LOG2RPKM=None,

                              to_smooth=False,
                              max_min_window=0,
                              compute_per_state_scores=False):

    foreground_scores = dict((k, {}) for k in [OVERALL_SCORE] + (states if compute_per_state_scores else []))

    datasets_A = sorted(group_A_segmentations)
    datasets_B = sorted(group_B_segmentations)

    for chrom in random_chunks:
        if emissions:
            chrom_segmentations_A = [chrom_segmentation_to_list(group_A_segmentations[d][chrom], BIN_SIZE) for d in datasets_A]
            chrom_segmentations_B = [chrom_segmentation_to_list(group_B_segmentations[d][chrom], BIN_SIZE) for d in datasets_B]
        elif signal_dir:
            chrom_signal_A = read_ChromHMM_signal(chrom, signal_dir, datasets_A, total_reads_per_mark, BIN_SIZE, LOG2RPKM)
            chrom_signal_B = read_ChromHMM_signal(chrom, signal_dir, datasets_B, total_reads_per_mark, BIN_SIZE, LOG2RPKM)

        for chunk_start, chunk_end in random_chunks[chrom]:

            if emissions:
                chunk_segmentations_A = dict((d, seg[chunk_start:chunk_end]) for d, seg in izip(datasets_A, chrom_segmentations_A))
                chunk_segmentations_B = dict((d, seg[chunk_start:chunk_end]) for d, seg in izip(datasets_B, chrom_segmentations_B))

                chunk_scores = get_overall_and_per_state_diff_score_emissions(chrom,
                                                                    chunk_segmentations_A,
                                                                    chunk_segmentations_B,
                                                                    emissions,
                                                                    BIN_SIZE,
                                                                    states,
                                                                    compute_per_state_scores,
                                                                    max_min_window,
                                                                    background_chunk=True)

            elif signal_dir:
                chunk_scores = compute_signal_diff(chrom_signal_A, chrom_signal_B, chunk_start, chunk_end, compute_per_state_scores=compute_per_state_scores)

            if to_smooth:
                chunk_scores = dict((score_type, smooth(chunk_scores[score_type])) for score_type in chunk_scores)

            for score_type in foreground_scores:
                for s in chunk_scores[score_type]:
                    abs_score = abs(s)
                    if abs_score not in foreground_scores[score_type]:
                        foreground_scores[score_type][abs_score] = 0

                    foreground_scores[score_type][abs_score] += 1

    for score_type in foreground_scores:
        total_regions = sum(foreground_scores[score_type].itervalues())
        total_so_far = 0

        for s in sorted(foreground_scores[score_type], reverse=True):
            total_so_far += foreground_scores[score_type][s]
            foreground_scores[score_type][s] = total_so_far / float(total_regions)

    return dict((score_type, sorted(foreground_scores[score_type].items())) for score_type in foreground_scores)

def compute_background_model(args,
                             fdr_threshold,
                             chromosomes,
                             group_A_segmentations,
                             group_B_segmentations,
                             states,
                             emissions,
                             signal_dir,
                             total_reads_per_mark,
                             BIN_SIZE,
                             LOG2RPKM):

    # pick the random chunks

    RANDOM_CHUNK_LENGTH = 1000000 / BIN_SIZE
    N_RANDOM_CHUNKS = 100
    MAX_N_SHUFFLED_SAMPLES = 100

    chrom_lengths = [group_A_segmentations.values()[0][chrom][-1][1] / BIN_SIZE for chrom in chromosomes]
    total_genome_length = sum(chrom_lengths)

    random_chunks = {}

    for _ in xrange(N_RANDOM_CHUNKS):
        chunk_start = None
        chrom = None

        while chunk_start is None:
            chunk_start = random.randint(0, total_genome_length - RANDOM_CHUNK_LENGTH - 1)

            offset = 0
            chrom_idx = 0
            while offset + chrom_lengths[chrom_idx] <= chunk_start:
                offset += chrom_lengths[chrom_idx]
                chrom_idx += 1

            chunk_start = chunk_start - offset
            chrom = chromosomes[chrom_idx]
            if chunk_start + RANDOM_CHUNK_LENGTH >= chrom_lengths[chrom_idx] or \
                    any(overlap(chunk_start, chunk_start + RANDOM_CHUNK_LENGTH, cs, ce) > 0
                        for cs, ce in random_chunks.get(chrom, [])):
                chunk_start = None

        if chrom not in random_chunks:
            random_chunks[chrom] = []

        random_chunks[chrom].append((chunk_start, chunk_start + RANDOM_CHUNK_LENGTH))

    random_chunks = dict((chrom, sorted(random_chunks[chrom])) for chrom in random_chunks)

    if args.chrom_hmm_binarized:
        from permute_marks_in_binarized_ChromHMM_baseline import compute_background_scores_by_shuffling_marks
        if emissions:
            background_scores = compute_background_scores_by_shuffling_marks(args.chrom_hmm_binarized,
                                                                         group_A_segmentations.keys() + group_B_segmentations.keys(),
                                                                         args.chrom_hmm_model_path,
                                                                         random_chunks,
                                                                         n_group_A=len(group_A_segmentations),
                                                                         BIN_SIZE=BIN_SIZE,
                                                                         LOG2RPKM=LOG2RPKM,
                                                                         states=states,
                                                                         n_perms=MAX_N_SHUFFLED_SAMPLES,
                                                                         compute_per_state_scores=args.per_state_scores,
                                                                         emissions=emissions,
                                                                         max_min_window=args.max_min_window,
                                                                         to_smooth=args.smooth
                                                                         )
    elif emissions:
        echo('Learning significance threshold for at p-value:', fdr_threshold)
        background_scores = compute_background_scores_by_shuffling_segmentations(group_A_segmentations,
                                                                                 group_B_segmentations,
                                                                                 random_chunks,
                                                                                 states,
                                                                                 emissions=emissions,
                                                                                 signal_dir=signal_dir,
                                                                                 total_reads_per_mark=total_reads_per_mark,
                                                                                 n_perm=MAX_N_SHUFFLED_SAMPLES,
                                                                                 compute_per_state_scores=args.per_state_scores,
                                                                                 to_smooth=args.smooth,
                                                                                 max_min_window=args.max_min_window,
                                                                                 n_threads=args.n_threads,
                                                                                 BIN_SIZE=BIN_SIZE,
                                                                                 LOG2RPKM=LOG2RPKM)
    elif signal_dir:
        from permute_marks_in_binarized_ChromHMM_baseline import compute_background_scores_by_shuffling_marks_signal
        background_scores = compute_background_scores_by_shuffling_marks_signal(group_A_segmentations.keys() + group_B_segmentations.keys(),
                                                                                 random_chunks,
                                                                                 n_group_A=len(group_A_segmentations),
                                                                                 n_perms=MAX_N_SHUFFLED_SAMPLES,
                                                                                 BIN_SIZE=BIN_SIZE,
                                                                                 LOG2RPKM=LOG2RPKM,
                                                                                 max_min_window=args.max_min_window,
                                                                                 signal_dir=signal_dir,
                                                                                 total_reads_per_mark=total_reads_per_mark,
                                                                                 to_smooth=args.smooth)

    foreground_scores = compute_foreground_scores(group_A_segmentations,
                                                  group_B_segmentations,
                                                  random_chunks,
                                                  states,
                                                  emissions=emissions,
                                                  signal_dir=signal_dir,
                                                  total_reads_per_mark=total_reads_per_mark,
                                                  to_smooth=args.smooth,
                                                  max_min_window=args.max_min_window,
                                                  compute_per_state_scores=args.per_state_scores,
                                                  BIN_SIZE=BIN_SIZE,
                                                  LOG2RPKM=LOG2RPKM)

    background_model = dict((k, {0.0: 1}) for k in [OVERALL_SCORE] + (states if compute_per_state_scores else []))
    background_threshold = dict((k, None) for k in [OVERALL_SCORE] + (states if compute_per_state_scores else []))

    # determine FDR threshold
    for score_type in sorted(background_model):

        all_scores = foreground_scores[score_type]

        best_score_cutoff = 0
        best_fdr = 1

        score_idx = 0
        bgr_idx = 0

        while score_idx < len(all_scores):

            # find how many background scores are greater than the current real score
            while bgr_idx < len(background_scores[score_type]) and background_scores[score_type][bgr_idx][0] < all_scores[score_idx][0]:
                bgr_idx += 1

            if bgr_idx == len(background_scores[score_type]):
                false_positives = 1
            else:
                false_positives = background_scores[score_type][bgr_idx][1]

            all_positives = all_scores[score_idx][1]

            current_fdr = false_positives / all_positives

            if current_fdr < best_fdr:
                best_fdr = current_fdr
                best_score_cutoff = all_scores[score_idx][0]
            else:
                current_fdr = best_fdr

            background_model[score_type][all_scores[score_idx][0]] = current_fdr

            if current_fdr <= fdr_threshold and background_threshold[score_type] is None:
                background_threshold[score_type] = all_scores[score_idx][0]

            score_idx += 1

        echo(score_type,
             'Score cutoff at FDR', fdr_threshold, ':', background_threshold[score_type],
             'bgr model size:', len(background_model[score_type]))

        if background_threshold[score_type] is None:
            echo('No significant changes were found at FDR of', fdr_threshold, '. Best FDR=', best_fdr, 'cutoff=', best_score_cutoff)
            background_threshold[score_type] = best_score_cutoff

    return dict((score_type, sorted(background_model[score_type].iteritems())) for score_type in background_model), background_threshold

def read_emissions(emissions_fname):
    metric = {}
    with open(emissions_fname) as in_f:
        header = in_f.readline().strip().split('\t')
        for l in in_f:
            buf = l.strip().split('\t')
            metric['E' + buf[0]] = dict((m, float(v)) for m, v in zip(header[1:], buf[1:]))
    return metric


def calculate_total_reads(celltypes, signal_dir):
    echo('Counting the reads')
    total_reads_per_mark = {}
    for fname in os.listdir(signal_dir):
        if not fname.endswith('txt.gz'):
            continue

        with open_file(os.path.join(signal_dir, fname)) as in_f:
            ct, chrom = in_f.readline().strip().split()
            if ct not in celltypes:# or chrom not in chromosomes:
                continue

            echo(fname)
            marks = in_f.readline().strip().split()

            if ct not in total_reads_per_mark:
                total_reads_per_mark[ct] = dict((m, 0) for m in marks)

            for i, l in enumerate(in_f):
                for m, v in izip(marks, map(int, l.strip().split())):
                    total_reads_per_mark[ct][m] += v
    print pprint.pformat(total_reads_per_mark)
    return total_reads_per_mark


def get_celltype(fname):
    return re.sub(r'_(\d+)(_coreMarks)?$', '', os.path.split(fname)[1].split('_segments')[0])


EMISSION_SCORES = 'emissions'
MARK_SCORES = 'marks'
POSTERIOR_SCORES = 'posteriors'
STATE_COUNT = 'state_count'
FISHER_EXACT = 'fisher'

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-a', nargs='+', help='segmentation files of group A')
    parser.add_argument('-b', nargs='+', help='segmentation files of group B')

    parser.add_argument('-A', dest='A', help='text file with segmentation files of group A, one per line')
    parser.add_argument('-B', dest='B', help='text file with segmentation files of group B, one per line')

    parser.add_argument('--bin-size', type=int, help='Bin size, default: %(default)s', default=200)

    parser.add_argument('-f', '--fdr-threshold', type=float, help='FDR threshold: %(default)s', default=0.01)

    parser.add_argument('-o', '--output', help='output prefix')

    parser.add_argument('--chrom-hmm-binarized', help='path to ChromHMM binarized files to compute background scores')
    parser.add_argument('--chrom-hmm-model-path', help='path to the ChromHMM model to use for shuffled marks')

    parser.add_argument('-s', '--smooth', action='store_true', help='smooth signal')
    parser.add_argument('-p', '--n-threads', type=int, help='number of threads to use', default=1)

    parser.add_argument('--max-min-window', type=int, help='pick the maximum distance between the two most similar bins within this window', default=0)

    parser.add_argument('--per-state-scores', action='store_true', help='compute per state scores', default=False)

    parser.add_argument('--emissions', help='emission probabilities to compute score based on euclidean distance of emissions')

    parser.add_argument('--signal-dir', help='ChromHMM signal directory for euclidean distance of the signal')
    parser.add_argument('--log2rpkm', action='store_true', help='take log2 rpkm for signal delta', default=False)
    parser.add_argument('--keep-max-scores', action='store_true', help='keep only the maximum score for each bin in each direction', default=False)
    parser.add_argument('--use_symmetric_KL', action='store_true', help='Use symmetric KL for posteriors delta, else Hellinger distance', default=False)
    parser.add_argument('--skip-summary', action='store_true', help='Skip summary bed files', default=False)

    parser.add_argument('--posteriors-dir', help='ChromHMM posteriors directory for computing the difference based on posteriors only')
    parser.add_argument('--baseline', choices=[EMISSION_SCORES, MARK_SCORES, POSTERIOR_SCORES, STATE_COUNT, FISHER_EXACT], help='baseline type')

    parser.add_argument('--filter-chromosomes', nargs='+', dest='filter_chroms',
                        help='Apply method only on these chromosomes. Default: use all chromosomes',
                        default=None)

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)

    open_log(args.output + '.log')
    BIN_SIZE = args.bin_size
    echo('cmd:', ' '.join(sys.argv))
    for arg in sorted(vars(args)):
        echo(arg, '=', getattr(args, arg))

    random.seed(42)

    if args.baseline is None:
        echo('--baseline is required parameter!')
        exit(1)

    LOG2RPKM = args.log2rpkm

    max_min_window = args.max_min_window

    compute_per_state_scores = args.per_state_scores
    keep_max_scores = args.keep_max_scores
    if args.baseline != FISHER_EXACT and keep_max_scores:
        error('--keep-max-scores is implemented only for Fisher exact test')

    if args.A:
        seg_A_fnames = [os.path.join(os.path.split(args.A)[0], l.strip()) for l in open_file(args.A)]
    else:
        seg_A_fnames = list(args.a)

    if args.B:
        seg_B_fnames = [os.path.join(os.path.split(args.B)[0], l.strip()) for l in open_file(args.B)]
    else:
        seg_B_fnames = list(args.b)

    group_A_segmentations, states = read_segmentations(seg_A_fnames)
    group_B_segmentations, _ = read_segmentations(seg_B_fnames)

    if args.filter_chroms is not None:
        chromosomes = args.filter_chroms
        echo('Using chromosomes:', chromosomes)
    else:
        chromosomes = sorted(reduce(operator.and_, [set(s.keys())
                                                 for segs in [group_A_segmentations, group_B_segmentations]
                                                    for s in segs.values()]),
                                    key=lambda c: int(c.replace('chr','')) if re.search(r'^chr(\d+)$', c) else 100)


    # chromosomes = sorted(reduce(operator.and_, [set(s.keys())
    #                                          for segs in [group_A_segmentations, group_B_segmentations]
    #                                          for s in segs.values()]),
    #                      key=lambda c: int(c.replace('chr','')) if re.search(r'^chr(\d+)$', c) else 100)
    # chromosomes = ['chr11']
    # print chromosomes
    group_A_segmentations = filter_chroms(group_A_segmentations, chromosomes)
    group_B_segmentations = filter_chroms(group_B_segmentations, chromosomes)

    emissions = None
    if args.emissions:
        emissions = read_emissions(args.emissions)

    posteriors_dir = None
    if args.posteriors_dir:
        posteriors_dir = args.posteriors_dir

    signal_dir = None
    total_reads_per_mark = None
    if args.signal_dir:
        signal_dir = args.signal_dir
        celltypes = sorted([get_celltype(d) for d in seg_A_fnames + seg_B_fnames])

        total_reads_per_mark = calculate_total_reads(celltypes, args.signal_dir)

    fdr_threshold = args.fdr_threshold

    out_fname = args.output

    # background_model, background_threshold = compute_background_model(args,
    #                                                                   fdr_threshold,
    #                                                                   chromosomes,
    #                                                                   group_A_segmentations,
    #                                                                   group_B_segmentations,
    #                                                                   states,
    #                                                                   emissions=emissions,
    #                                                                   signal_dir=signal_dir,
    #                                                                   total_reads_per_mark=total_reads_per_mark,
    #                                                                   BIN_SIZE=BIN_SIZE,
    #                                                                   LOG2RPKM=LOG2RPKM)

    real_scores = {}
    real_state_transitions = {}
    echo('Computing real scores')

    # pure_state_transitions = generate_pure_state_transitions(emissions=emissions)
    # chromosomes = ['chrX']
    if args.n_threads > 1:
        pool = Pool(args.n_threads)
        _map = pool.imap_unordered

    all_signal = {}
    all_transitions = {}
    score_types = None
    for chrom, signal, state_transitions in _map(worker,
                                                 [{EMISSION_SCORES:
                                                       (get_overall_and_per_state_diff_score_emissions,
                                                       _chrom,
                                                       dict((d, {_chrom: group_A_segmentations[d][_chrom]}) for d in group_A_segmentations),
                                                       dict((d, {_chrom: group_B_segmentations[d][_chrom]}) for d in group_B_segmentations),
                                                       emissions if emissions else signal_dir,
                                                       BIN_SIZE,
                                                       total_reads_per_mark,
                                                       compute_per_state_scores,
                                                       max_min_window),
                                                   MARK_SCORES:(get_overall_and_per_state_diff_score_signal,
                                                                _chrom,
                                                                dict((d, {_chrom: group_A_segmentations[d][_chrom]}) for d in group_A_segmentations),
                                                                dict((d, {_chrom: group_B_segmentations[d][_chrom]}) for d in group_B_segmentations),
                                                                signal_dir,
                                                                BIN_SIZE,
                                                                LOG2RPKM,
                                                                total_reads_per_mark,
                                                                compute_per_state_scores,
                                                                max_min_window),
                                                   POSTERIOR_SCORES:(
                                                        get_overall_and_per_state_diff_score_posteriors,
                                                        _chrom,
                                                        dict((d, {_chrom: group_A_segmentations[d][_chrom]}) for d in group_A_segmentations),
                                                        dict((d, {_chrom: group_B_segmentations[d][_chrom]}) for d in group_B_segmentations),
                                                        posteriors_dir,
                                                        BIN_SIZE,
                                                        compute_per_state_scores,
                                                        args.use_symmetric_KL
                                                   ),
                                                   STATE_COUNT:(
                                                        get_overall_and_per_state_diff_score_by_state_counts,
                                                        _chrom,
                                                        dict((d, {_chrom: group_A_segmentations[d][_chrom]}) for d in group_A_segmentations),
                                                        dict((d, {_chrom: group_B_segmentations[d][_chrom]}) for d in group_B_segmentations),
                                                        BIN_SIZE,
                                                        states,
                                                        compute_per_state_scores
                                                   ),
                                                   FISHER_EXACT:(
                                                        get_overall_and_per_state_diff_score_by_fishers_exact_test,
                                                        _chrom,
                                                        dict((d, {_chrom: group_A_segmentations[d][_chrom]}) for d in group_A_segmentations),
                                                        dict((d, {_chrom: group_B_segmentations[d][_chrom]}) for d in group_B_segmentations),
                                                        BIN_SIZE,
                                                        states,
                                                        compute_per_state_scores,
                                                        keep_max_scores
                                                   ),
                                                   }[args.baseline] for _chrom in chromosomes]):

        if args.smooth:
            smooth_dict(signal)

        score_types = sorted(signal)
        print score_types
        all_signal[chrom] = signal
        all_transitions[chrom] = state_transitions

    if args.n_threads > 1:
        pool.close()

    for score_type in score_types:
        gc.collect()
        echo('Writing wiggle file for', score_type)
        with open_file(out_fname + '.' + re.sub(r'\W+', '_', score_type) + '.wig.gz', 'w') as out_signal_file:

            title = os.path.split(out_fname)[1].replace('.wig', '').replace('.gz', '') + ' ' + score_type
            out_signal_file.write('track type=wiggle_0 name="%s" description="%s"\n' % (title, title))

            for chrom in all_signal:
                store_wig(chrom, all_signal[chrom][score_type], out_signal_file, BIN_SIZE)

        if not args.skip_summary:
            echo('Writing sorted summary files for', score_type)

            def write_line(out_f, chrom, bin_no, score, closest_state_A, closest_state_B, a_states, b_states, score_type, threshold, adjusted_pvalue='NA'):
                out_f.write('\t'.join(map(str,
                                          [ chrom,
                                            bin_no * BIN_SIZE,
                                            (bin_no + 1) * BIN_SIZE,
                                            '.',
                                            score,
                                            '.',
                                            closest_state_A + '-' + closest_state_B,
                                            ','.join(a_states) + '/' + ','.join(b_states),
                                            adjusted_pvalue, #get_FDR(abs(score), background_model, score_type),
                                            '1' if score >= threshold else '0'])) + '\n')


            group_scores = [(score, chrom, bin_no, state_transition)
                            for chrom in chromosomes
                            for bin_no, (score, state_transition) in
                            enumerate(izip(all_signal[chrom][score_type], all_transitions[chrom]))]

            if score_type is OVERALL_SCORE:
                sorted_scores = sorted(group_scores,
                                       key=lambda (score, chrom, bin_no, state_transition): (-score, chrom, bin_no))

                bed_out_fname = out_fname + '.' + re.sub(r'\W+', '_', score_type) + '.summary.bed.gz'

            else:
                sorted_scores = sorted(group_scores,
                                      key=lambda (score, chrom, bin_no, state_transition): (score, chrom, bin_no))

                bed_out_fname = out_fname + '.' + re.sub(r'\W+', '_', score_type) + '.both_groups.summary.bed.gz'

            if args.baseline == FISHER_EXACT:
                from statsmodels.sandbox.stats.multicomp import multipletests
                echo('Computing adjusted p-values')
                pvals = [2**(-abs(s)) for s, _, _, _ in sorted_scores]
                _, adjusted_pvals, _, _ = multipletests(pvals, method='fdr_bh')
            else:
                adjusted_pvals = None

            with open_file(bed_out_fname, 'w') as out_f:
                for score_idx, (score, chrom, bin_no, (closest_state_A, closest_state_B, a_states, b_states)) in enumerate(sorted_scores):
                    write_line(out_f,
                               chrom,
                               bin_no,
                               score,
                               closest_state_A,
                               closest_state_B,
                               a_states,
                               b_states,
                               score_type,
                               0,
                               adjusted_pvalue=adjusted_pvals[score_idx] if adjusted_pvals is not None else 'NA')
        #
        # if score_type is OVERALL_SCORE:
        #     group_scores = [(score, chrom, bin_no, state_transition)
        #                     for chrom in chromosomes
        #                         for bin_no, (score, state_transition) in enumerate(izip(all_signal[chrom][score_type], all_transitions[chrom]))]
        #
        #     with open_file(out_fname + '.' + re.sub(r'\W+', '_', score_type) + '.summary.bed.gz', 'w') as out_f:
        #         for score, chrom, bin_no, (closest_state_A, closest_state_B, a_states, b_states) \
        #                 in sorted(group_scores,
        #                           key=lambda (score, chrom, bin_no, state_transition): (-score, chrom, bin_no)):
        #
        #             write_line(out_f, chrom, bin_no, score, closest_state_A, closest_state_B, a_states, b_states, score_type, 0) #background_threshold[score_type])
        #
        # else:
        #     group_scores = [(score, chrom, bin_no, state_transition)
        #                     for chrom in chromosomes
        #                         for bin_no, (score, state_transition) in enumerate(izip(all_signal[chrom][score_type], all_transitions[chrom]))]
        #
        #     with open_file(out_fname + '.' + re.sub(r'\W+', '_', score_type) + '.both_groups.summary.bed.gz', 'w') as out_f:
        #         for score, chrom, bin_no, (closest_state_A, closest_state_B, a_states, b_states) \
        #                 in sorted(group_scores,
        #                           key=lambda (score, chrom, bin_no, state_transition): (score, chrom, bin_no)):
        #
        #             write_line(out_f, chrom, bin_no, score, closest_state_A, closest_state_B, a_states, b_states, score_type, 0) #background_threshold[score_type])

            # for group in [GROUP_A, GROUP_B]:
            #     group_scores = [(score, chrom, bin_no, state_transition)
            #                     for chrom in chromosomes
            #                         for bin_no, (score, state_transition) in enumerate(izip(all_signal[chrom][score_type], all_transitions[chrom]))
            #                             if ((group == GROUP_A and score >= 0) or (group == GROUP_B and score < 0))]
            #
            #     with open_file(out_fname + '.' + re.sub(r'\W+', '_', score_type) + '.group_' + group + '.summary.bed.gz', 'w') as out_f:
            #         for score, chrom, bin_no, (closest_state_A, closest_state_B, a_states, b_states) \
            #                 in sorted(group_scores,
            #                           key=lambda (score, chrom, bin_no, state_transition): (-abs(score), chrom, bin_no)):
            #             write_line(out_f, chrom, bin_no, score, closest_state_A, closest_state_B, a_states, b_states, score_type, 0) #background_threshold[score_type])
            #
        for chrom in all_signal:
            del all_signal[chrom][score_type]

    echo('Done')
    close_log()
