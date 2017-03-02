import argparse
from itertools import izip, combinations
from multiprocessing import Pool
import os
import random
import sys
import operator
import pickle
import itertools
import gc
import traceback

from utils import open_file, echo, n_combinations, smooth
from metric import *



def compress_segmentations_into_non_repeating_windows(chrom_segmentations_A, chrom_segmentations_B, max_min_window):

    chrom_length = len(chrom_segmentations_A[0])

    non_repating_windows = []
    prev_a_windows = None
    prev_b_windows = None
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

        if a_window != prev_a_windows or b_window != prev_b_windows:

            key = (a_bins, b_bins, a_window, b_window)

            if key not in unique_windows:
                unique_windows_indices[key] = len(unique_windows)
                unique_windows.add(key)

            non_repating_windows.append([unique_windows_indices[key], 0])

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

def get_overall_and_per_state_diff_score(chrom,
                                         group_A_segmentations,
                                         metric_A,
                                         group_B_segmentations,
                                         metric_B,
                                         BIN_SIZE,
                                         states,
                                         pure_state_transitions,
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
    gc.collect()

    if not background_chunk:
        echo('Done compressing into windows:', len(chrom_segmentations_windows), 'unique windows:', len(unique_windows))

    chrom_length = len(chrom_segmentations_windows)
    scores = [None] * chrom_length # dict((k, [0.0] * chrom_length) for k in [OVERALL_SCORE] + (states if compute_per_state_scores else []))

    scores_cache = [None] * len(unique_windows)
    closest_transition_cache = [None] * len(unique_windows)

    state_transitions = None
    if not background_chunk:
        state_transitions = [None] * chrom_length

    # return scores, state_transitions

    for bin_idx, (window_idx, _) in enumerate(chrom_segmentations_windows):

        a_bins, b_bins, a_window, b_window = unique_windows[window_idx]

        if scores_cache[window_idx] is None:
            cache_miss += 1

            # get the posterior probabilities for group A by summing up the
            # the probability vectors for each state
            a_probs = get_prob_vector(a_bins, datasets_A, metric_A, states)

            # same for group B
            b_probs = get_prob_vector(b_bins, datasets_B, metric_B, states)

            a_window_probs = [get_prob_vector(a_window[i * len(datasets_A):(i + 1) * len(datasets_A)], datasets_A, metric_A, states)
                              for i in xrange(len(a_window) / len(datasets_A))]

            b_window_probs = [get_prob_vector(b_window[i * len(datasets_B):(i + 1) * len(datasets_B)], datasets_B, metric_B, states)
                              for i in xrange(len(b_window) / len(datasets_B))]

            score, (true_a_probs, true_b_probs, true_a_idx, true_b_idx) = max(
                min([(symmetric_KL_divergence(a_probs, probs), (a_probs, probs, -1, b_probs_i))
                     for b_probs_i, probs in enumerate(b_window_probs)]),
                min([(symmetric_KL_divergence(probs, b_probs), (probs, b_probs, a_probs_i, -1))
                     for a_probs_i, probs in enumerate(a_window_probs)]))

            scores_dict = dict((k, 0) for k in [OVERALL_SCORE] + (states if compute_per_state_scores else []))
            scores_dict[OVERALL_SCORE] = score

            if compute_per_state_scores:
                for state in states:
                    scores_dict[state] = true_a_probs[state] - true_b_probs[state]

            if not background_chunk:
                flat_prob_dict = dict((group + ' ' + key, d[key])
                                      for group, d in [('A', true_a_probs), ('B', true_b_probs)] for key in d)

                true_a_bins = set(a_bins if true_a_idx == -1 else a_window[true_a_idx * len(datasets_A):(true_a_idx + 1) * len(datasets_A)])
                true_b_bins = set(b_bins if true_b_idx == -1 else b_window[true_b_idx * len(datasets_B):(true_b_idx + 1) * len(datasets_B)])

                state_pairs = [(s1, s2) for s1 in true_a_bins for s2 in true_b_bins]

                closest_state_A, closest_state_B = min(state_pairs,
                                                       key=lambda (_s1, _s2):
                                                       KL(flat_prob_dict,
                                                          pure_state_transitions[_s1][_s2][STATE_DISTRIBUTION]))

                closest_transition_cache[window_idx] = (closest_state_A, closest_state_B, a_bins, b_bins)

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

    for bin_no, bin_signal in enumerate(signal):

        out_signal_f.write('%.2lf\n' % bin_signal)


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


def _compute_background_scores_by_shuffling_segmentations(group_A_segmentations,
                                                          group_B_segmentations,
                                                          random_chunks,
                                                          states,
                                                          metric_A,
                                                          metric_B,
                                                          n_perm=100,
                                                          to_smooth=False,
                                                          max_min_window=0,
                                                          compute_per_state_scores=False
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

    combo_no = 0

    background_scores = dict((k, {}) for k in [OVERALL_SCORE] + (states if compute_per_state_scores else []))

    for shuffled_A_datasets, shuffled_B_datasets in shuffled_samples(n_perm):

        combo_no += 1
        echo('Combo no:', combo_no)

        print 'Shuffled A:', shuffled_A_datasets
        print 'Shuffled B:', shuffled_B_datasets
        print

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

        # iterate over the bins in both segmentation groups
        # for chrom in random.sample(chromosomes, 2):

        for chrom in random_chunks:

            chrom_segmentations_A = [chrom_segmentation_to_list(shuffled_A_segmentations[d][chrom], BIN_SIZE) for d in shuffled_A_datasets]
            chrom_segmentations_B = [chrom_segmentation_to_list(shuffled_B_segmentations[d][chrom], BIN_SIZE) for d in shuffled_B_datasets]

            for chunk_start, chunk_end in random_chunks[chrom]:
                chunk_segmentations_A = dict((d, seg[chunk_start:chunk_end]) for d, seg in izip(shuffled_A_segmentations, chrom_segmentations_A))
                chunk_segmentations_B = dict((d, seg[chunk_start:chunk_end]) for d, seg in izip(shuffled_B_segmentations, chrom_segmentations_B))

                chunk_scores = get_overall_and_per_state_diff_score(chrom,
                                                                    chunk_segmentations_A,
                                                                    shuffled_metric_A,
                                                                    chunk_segmentations_B,
                                                                    shuffled_metric_B,
                                                                    BIN_SIZE,
                                                                    states,
                                                                    None,
                                                                    compute_per_state_scores,
                                                                    max_min_window,
                                                                    background_chunk=True)

                if to_smooth:
                    chunk_scores = dict((score_type, smooth(chunk_scores[score_type])) for score_type in chunk_scores)

                for score_type in background_scores:
                    for s in chunk_scores[score_type]:
                        abs_score = abs(s)
                        if abs_score not in background_scores[score_type]:
                            background_scores[score_type][abs_score] = 0

                        background_scores[score_type][abs_score] += 1

    for score_type in background_scores:
        total_regions = sum(background_scores[score_type].itervalues())
        total_so_far = 0

        for s in sorted(background_scores[score_type], reverse=True):
            total_so_far += background_scores[score_type][s]
            background_scores[score_type][s] = total_so_far / float(total_regions)

    return dict((score_type, sorted(background_scores[score_type].items())) for score_type in background_scores)


def compute_background_scores_by_shuffling_segmentations(group_A_segmentations,
                                                         group_B_segmentations,
                                                         random_chunks,
                                                         states,
                                                         metric_A,
                                                         metric_B,
                                                         n_perm=100,
                                                         to_smooth=False,
                                                         max_min_window=0,
                                                         compute_per_state_scores=False,
                                                         n_threads=1
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
        _map = pool.map

    for combo_background_scores in _map(worker,
                                        [(process_shuffled_combo,
                                          combo_no,
                                          shuffled_A_datasets,
                                          dict((dataset, all_segs[dataset]) for dataset in shuffled_A_datasets),
                                          dict((dataset, all_metrics[dataset]) for dataset in shuffled_A_datasets),

                                          shuffled_B_datasets,
                                          dict((dataset, all_segs[dataset]) for dataset in shuffled_B_datasets),
                                          dict((dataset, all_metrics[dataset]) for dataset in shuffled_B_datasets),

                                          random_chunks,
                                          BIN_SIZE,
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


def process_shuffled_combo( combo_no,
                            shuffled_A_datasets,
                            shuffled_A_segmentations,
                            shuffled_metric_A,

                            shuffled_B_datasets,
                            shuffled_B_segmentations,
                            shuffled_metric_B,

                            random_chunks,
                            BIN_SIZE,
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
        chrom_segmentations_A = [chrom_segmentation_to_list(shuffled_A_segmentations[d][chrom], BIN_SIZE) for d in shuffled_A_datasets]
        chrom_segmentations_B = [chrom_segmentation_to_list(shuffled_B_segmentations[d][chrom], BIN_SIZE) for d in shuffled_B_datasets]

        for chunk_start, chunk_end in random_chunks[chrom]:
            # echo(combo_no, chunk_start, chunk_end)
            chunk_segmentations_A = dict((d, seg[chunk_start:chunk_end]) for d, seg in izip(shuffled_A_segmentations, chrom_segmentations_A))
            chunk_segmentations_B = dict((d, seg[chunk_start:chunk_end]) for d, seg in izip(shuffled_B_segmentations, chrom_segmentations_B))

            chunk_scores = get_overall_and_per_state_diff_score(chrom,
                                                                chunk_segmentations_A,
                                                                shuffled_metric_A,
                                                                chunk_segmentations_B,
                                                                shuffled_metric_B,
                                                                BIN_SIZE,
                                                                states,
                                                                None,
                                                                compute_per_state_scores,
                                                                max_min_window,
                                                                background_chunk=True)

            if to_smooth:
                chunk_scores = dict((score_type, smooth(chunk_scores[score_type])) for score_type in chunk_scores)

            for score_type in background_scores:
                for s in chunk_scores[score_type]:
                    abs_score = abs(s)
                    if abs_score not in background_scores[score_type]:
                        background_scores[score_type][abs_score] = 0

                    background_scores[score_type][abs_score] += 1

    return background_scores

def compute_foreground_scores(group_A_segmentations,
                              group_B_segmentations,
                              random_chunks,
                              states,
                              metric_A,
                              metric_B,
                              to_smooth=False,
                              max_min_window=0,
                              compute_per_state_scores=False):

    foreground_scores = dict((k, {}) for k in [OVERALL_SCORE] + (states if compute_per_state_scores else []))

    datasets_A = sorted(group_A_segmentations)
    datasets_B = sorted(group_B_segmentations)

    for chrom in random_chunks:

        chrom_segmentations_A = [chrom_segmentation_to_list(group_A_segmentations[d][chrom], BIN_SIZE) for d in datasets_A]
        chrom_segmentations_B = [chrom_segmentation_to_list(group_B_segmentations[d][chrom], BIN_SIZE) for d in datasets_B]

        for chunk_start, chunk_end in random_chunks[chrom]:
            chunk_segmentations_A = dict((d, seg[chunk_start:chunk_end]) for d, seg in izip(datasets_A, chrom_segmentations_A))
            chunk_segmentations_B = dict((d, seg[chunk_start:chunk_end]) for d, seg in izip(datasets_B, chrom_segmentations_B))

            chunk_scores = get_overall_and_per_state_diff_score(chrom,
                                                                chunk_segmentations_A,
                                                                metric_A,
                                                                chunk_segmentations_B,
                                                                metric_B,
                                                                BIN_SIZE,
                                                                states,
                                                                None,
                                                                compute_per_state_scores,
                                                                max_min_window,
                                                                background_chunk=True)

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
                             metric_A,
                             group_B_segmentations,
                             metric_B,
                             states):

    # pick the random chunks
    BIN_SIZE = args.bin_size

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
        from permute_marks_in_binarized_ChromHMM import compute_background_scores_by_shuffling_marks
        background_scores = compute_background_scores_by_shuffling_marks(args.chrom_hmm_binarized,
                                                                         group_A_segmentations.keys() + group_B_segmentations.keys(),
                                                                         args.chrom_hmm_model_path,
                                                                         random_chunks,
                                                                         n_group_A=len(group_A_segmentations),
                                                                         BIN_SIZE=BIN_SIZE,
                                                                         states=states,
                                                                         n_perms=MAX_N_SHUFFLED_SAMPLES,
                                                                         compute_per_state_scores=args.per_state_scores,
                                                                         max_min_window=args.max_min_window,
                                                                         to_smooth=args.smooth)

    else:
        echo('Learning significance threshold for at p-value:', fdr_threshold)
        background_scores = compute_background_scores_by_shuffling_segmentations(group_A_segmentations,
                                                                                 group_B_segmentations,
                                                                                 random_chunks,
                                                                                 states,
                                                                                 metric_A,
                                                                                 metric_B,
                                                                                 n_perm=MAX_N_SHUFFLED_SAMPLES,
                                                                                 compute_per_state_scores=args.per_state_scores,
                                                                                 to_smooth=args.smooth,
                                                                                 max_min_window=args.max_min_window,
                                                                                 n_threads=args.n_threads)

    foreground_scores = compute_foreground_scores(group_A_segmentations,
                                                  group_B_segmentations,
                                                  random_chunks,
                                                  states,
                                                  metric_A,
                                                  metric_B,
                                                  to_smooth=False,
                                                  max_min_window=0,
                                                  compute_per_state_scores=args.per_state_scores)

    background_model = dict((k, {0.0: 1}) for k in [OVERALL_SCORE] + (states if compute_per_state_scores else []))
    background_threshold = dict((k, None) for k in [OVERALL_SCORE] + (states if compute_per_state_scores else []))

    # determine FDR threshold
    for score_type in background_model:

        all_scores = foreground_scores[score_type]

        best_score_cutoff = 0
        best_fdr = 1

        score_idx = 0
        bgr_idx = 0

        while True:

            # if we reached the end of all_scores, then nothing is significant
            if score_idx == len(all_scores):
                break

            # find how many background scores are greater than the current real score
            while bgr_idx < len(background_scores[score_type]) and background_scores[score_type][bgr_idx][0] < all_scores[score_idx][0]:
                bgr_idx += 1

            if bgr_idx == len(background_scores[score_type]):
                false_positives = 0
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

        echo(score_type, 'Score cutoff at FDR', fdr_threshold, ':', background_threshold[score_type])

        if background_threshold[score_type] is None:
            echo('No significant changes were found at FDR of', fdr_threshold, '. Best FDR=', best_fdr, 'cutoff=', best_score_cutoff)
            background_threshold[score_type] = best_score_cutoff

    return dict((score_type, sorted(background_model[score_type].iteritems())) for score_type in background_model), background_threshold



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-a', nargs='+', help='segmentation files of group A')
    parser.add_argument('-b', nargs='+', help='segmentation files of group B')
    parser.add_argument('--bin-size', type=int, help='Bin size, default: %(default)s', default=200)
    parser.add_argument('-f', '--fdr-threshold', type=float, help='FDR threshold: %(default)s', default=0.01)

    parser.add_argument('-o', '--output', help='output prefix')
    parser.add_argument('--background-scores', help='bacground_scores.pickle')

    parser.add_argument('--chrom-hmm-binarized', help='path to ChromHMM binarized files to compute background scores')
    parser.add_argument('--chrom-hmm-model-path', help='path to the ChromHMM model to use for shuffled marks')

    parser.add_argument('-s', '--smooth', action='store_true', help='smooth signal')
    parser.add_argument('-p', '--n-threads', type=int, help='number of threads to use', default=1)

    parser.add_argument('--max-min-window', type=int, help='pick the maximum distance between the two most similar bins within this window', default=0)

    parser.add_argument('--per-state-scores', action='store_true', help='compute per state scores', default=False)
    # parser.add_argument('--use-closest-rep', action='store_true', help='use only closest replicate to learn the metric')
    # parser.add_argument('--all-for-null', action='store_true', default=False, help='include the original combo in the background model')
    parser.add_argument('--store-scores', action='store_true', default=False, help='store foreground and background scores in a pickle file')
    # parser.add_argument('--output-changes', action='store_true', default=False, help='output significant pairwise changes')

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)

    open_log(args.output + '.log')
    BIN_SIZE = args.bin_size
    echo('cmd:', ' '.join(sys.argv))
    for arg in sorted(vars(args)):
        echo(arg , '=', getattr(args, arg))

    random.seed(36830804669286)

    max_min_window = args.max_min_window

    compute_per_state_scores = args.per_state_scores

    group_A_segmentations, states = read_segmentations(args.a)
    group_B_segmentations, _ = read_segmentations(args.b)

    chromosomes = sorted(reduce(operator.and_, [set(s.keys())
                                             for segs in [group_A_segmentations, group_B_segmentations]
                                             for s in segs.values()]),
                         key=lambda c: int(c.replace('chr','')) if re.search(r'^chr(\d+)$', c) else 100)
    # print chromosomes
    group_A_segmentations = filter_chroms(group_A_segmentations, chromosomes)
    group_B_segmentations = filter_chroms(group_B_segmentations, chromosomes)

    metric_A = learn_metric_from_all_replicates(group_A_segmentations, states, BIN_SIZE, other_group=group_B_segmentations)
    print_metric(metric_A)

    metric_B = learn_metric_from_all_replicates(group_B_segmentations, states, BIN_SIZE, other_group=group_A_segmentations)
    print_metric(metric_B)
    print_average_metric(states, metric_A, metric_B)

    fdr_threshold = args.fdr_threshold

    out_fname = args.output

    background_model, background_threshold = compute_background_model(args,
                                                                      fdr_threshold,
                                                                      chromosomes,
                                                                      group_A_segmentations,
                                                                      metric_A,
                                                                      group_B_segmentations,
                                                                      metric_B,
                                                                      states)

    real_scores = {}
    real_state_transitions = {}
    echo('Computing real scores')

    out_signal_files = dict((score_type, open_file(out_fname + '.' + re.sub(r'\W+', '_', score_type) + '.wig.gz', 'w'))
                            for score_type in [OVERALL_SCORE] + (states if compute_per_state_scores else []))

    out_summary_files = dict((score_type, open_file(out_fname + '.' + re.sub(r'\W+', '_', score_type) + '.summary.bed.gz', 'w'))
                            for score_type in [OVERALL_SCORE] + (states if compute_per_state_scores else []))

    for score_type in out_signal_files:
        title = os.path.split(out_fname)[1].replace('.wig', '').replace('.gz', '') + ' ' + score_type
        out_signal_files[score_type].write('track type=wiggle_0 name="%s" description="%s"\n' % (title, title))

    pure_state_transitions = generate_pure_state_transitions(states, metric_A, metric_B)
    # chromosomes = ['chrX']
    if args.n_threads > 1:
        pool = Pool(args.n_threads)
        _map = pool.map
    else:
        _map = itertools.imap

    for chrom, signal, state_transitions in _map(worker,
                                                 [(get_overall_and_per_state_diff_score,
                                                   _chrom,
                                                   dict((d, {_chrom: group_A_segmentations[d][_chrom]}) for d in group_A_segmentations),
                                                   metric_A,
                                                   dict((d, {_chrom: group_B_segmentations[d][_chrom]}) for d in group_B_segmentations),
                                                   metric_B,
                                                   BIN_SIZE,
                                                   states,
                                                   pure_state_transitions,
                                                   compute_per_state_scores,
                                                   max_min_window) for _chrom in chromosomes]):

    # for chrom in chromosomes:
    #     signal, state_transitions = get_overall_and_per_state_diff_score(chrom,
    #                                                                      group_A_segmentations,
    #                                                                      metric_A,
    #                                                                      group_B_segmentations,
    #                                                                      metric_B,
    #                                                                      BIN_SIZE,
    #                                                                      states,
    #                                                                      pure_state_transitions,
    #                                                                      compute_per_state_scores,
    #                                                                      max_min_window)

        if args.smooth:
            smooth_dict(signal)

        echo('Storing output files')

        store_dict_bed(chrom,
                       state_transitions,
                       signal,
                       background_model,
                       out_summary_files,
                       span=BIN_SIZE)

        echo('Storing wiggle')
        store_dict_wig(chrom, signal, out_signal_files, span=BIN_SIZE)

        gc.collect()

    if args.n_threads > 1:
        pool.close()

    for score_type in out_signal_files:
        out_signal_files[score_type].close()


    # with open_file(out_fname + '.summary.bed.gz', 'w') as out_states_f:
    #     for chrom in chromosomes:
    #         store_bed(chrom,
    #                   real_state_transitions[chrom],
    #                   real_scores[chrom],
    #                   score_to_fdr,
    #                   states,
    #                   out_states_f,
    #                   compute_per_state_scores)
    #
    #
    #
    # else:
    #     with open_file(out_fname + '.differential_peaks_FDR_' + str(fdr_threshold) + '.bed.gz', 'w') as out_diff_peaks_f:
    #         peak_no = 0
    #
    #         for chrom in sorted(real_scores):
    #
    #             # peak_start = None
    #             #
    #             # peak_score = 0
    #             # peak_transition = None
    #
    #             for bin_no, (bin_signal, (closest_state_A, closest_state_B, a_states, b_states)) in \
    #                     enumerate(izip(real_scores[chrom][OVERALL_SCORE],
    #                                    real_state_transitions[chrom])):
    #
    #                 bin_start = bin_no * BIN_SIZE
    #
    #                 if bin_signal >= score_cutoff:
    #                     peak_no += 1
    #                     out_diff_peaks_f.write('\t'.join(map(str, [chrom,
    #                                                                bin_start,
    #                                                                bin_start + BIN_SIZE,
    #                                                                'bin_' + str(peak_no),
    #                                                                '%.2lf' % bin_signal,
    #                                                                '.',
    #                                                                closest_state_A + '-' + closest_state_B])) + '\n')
    #
    #                 # if bin_signal >= score_cutoff:
    #                 #     if peak_start is None:
    #                 #         peak_start = bin_start
    #                 #         peak_score = 0
    #                 #         peak_no += 1
    #                 #
    #                 #     peak_score += bin_signal
    #                 #
    #                 # if (bin_signal < score_cutoff or bin_no == len(real_scores[chrom]) - 1) and peak_start is not None:
    #                 #
    #                 #     if bin_no == len(real_scores[chrom]) - 1:
    #                 #         bin_start += BIN_SIZE
    #                 #
    #                 #     out_diff_peaks_f.write('\t'.join(map(str, [chrom,
    #                 #                                                peak_start,
    #                 #                                                bin_start,
    #                 #                                                chrom + '_' + str(peak_no),
    #                 #                                                '%.2lf' % (BIN_SIZE * peak_score / (bin_start - peak_start)),
    #                 #                                                '.'])) + '\n')
    #                 #     peak_start = None
    #
    # if args.store_scores:
    #     echo('Storing scores in:', out_fname + '.scores.pickle')
    #     with open(out_fname + '.scores.pickle', 'w') as scores_f:
    #         pickle.dump({'bgr': background_scores,
    #                      'fgr': all_scores}, scores_f, protocol=pickle.HIGHEST_PROTOCOL)


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

    echo('Done')
    close_log()
