import argparse
from itertools import izip, combinations
from multiprocessing import Pool
import os
import random
import sys
import operator
import cPickle as pickle
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
    # echo('Expanding')
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


def get_overall_and_per_state_diff_score(chrom,
                                         group_A_segmentations,
                                         metric_A,
                                         group_B_segmentations,
                                         metric_B,
                                         states,
                                         consistent_state_probs,
                                         compute_per_state_scores,
                                         max_min_window,
                                         background_chunk=False,
                                         use_posteriors=False,
                                         report_per_group_state_probabilities=False,
                                         keep_max_scores_per_bin=False):

    verbose = not background_chunk
    # verbose = True

    if verbose:
        echo('Chromosome:', chrom)

    cache_miss = 0

    datasets_A = sorted(group_A_segmentations)
    datasets_B = sorted(group_B_segmentations)

    chrom_segmentations_A = [group_A_segmentations[d] for d in datasets_A]
    chrom_segmentations_B = [group_B_segmentations[d] for d in datasets_B]

    if verbose:
        echo('Chrom length:', len(chrom_segmentations_A[0]))

    chrom_segmentations_windows, unique_windows = compress_segmentations_into_non_repeating_windows(
                                                                                    chrom_segmentations_A,
                                                                                    chrom_segmentations_B,
                                                                                    max_min_window)

    del chrom_segmentations_A
    del chrom_segmentations_B

    if verbose:
        gc.collect()
        echo('Done compressing into windows:', len(chrom_segmentations_windows), 'unique windows:', len(unique_windows))

    chrom_length = len(chrom_segmentations_windows)
    scores = [None] * chrom_length # dict((k, [0.0] * chrom_length) for k in [OVERALL_SCORE] + (states if compute_per_state_scores else []))

    scores_cache = [None] * len(unique_windows)

    state_transitions = None
    closest_transition_cache = None

    if not background_chunk:
        state_transitions = [None] * chrom_length
        closest_transition_cache = [None] * len(unique_windows)

    # return scores, state_transitions
    if use_posteriors:
        _get_prob_vector = get_prob_vector_based_on_posteriors
    else:
        _get_prob_vector = get_prob_vector

    for bin_idx, (window_idx, _) in enumerate(chrom_segmentations_windows):

        a_bins, b_bins, a_window, b_window = unique_windows[window_idx]

        if scores_cache[window_idx] is None:
            cache_miss += 1

            # get the posterior probabilities for group A by summing up the
            # the probability vectors for each state
            a_probs = _get_prob_vector(a_bins, datasets_A, metric_A, states)

            # same for group B
            b_probs = _get_prob_vector(b_bins, datasets_B, metric_B, states)

            a_window_probs = [_get_prob_vector(a_window[i * len(datasets_A):(i + 1) * len(datasets_A)], datasets_A, metric_A, states)
                              for i in xrange(len(a_window) / len(datasets_A))]

            b_window_probs = [_get_prob_vector(b_window[i * len(datasets_B):(i + 1) * len(datasets_B)], datasets_B, metric_B, states)
                              for i in xrange(len(b_window) / len(datasets_B))]

            score, (true_a_probs, true_b_probs, true_a_idx, true_b_idx) = max(
                min([(symmetric_KL_divergence(a_probs, probs), (a_probs, probs, -1, b_probs_i))
                     for b_probs_i, probs in enumerate(b_window_probs)]),
                min([(symmetric_KL_divergence(probs, b_probs), (probs, b_probs, a_probs_i, -1))
                     for a_probs_i, probs in enumerate(a_window_probs)]))

            scores_dict = dict((k, 0) for k in [OVERALL_SCORE] + (states if compute_per_state_scores else []))
            scores_dict[OVERALL_SCORE] = round(score, PRECISION)

            if compute_per_state_scores:
                for state in states:
                    scores_dict[state] = round(true_a_probs[state] - true_b_probs[state], PRECISION)

                    if report_per_group_state_probabilities:
                        scores_dict[args.group_A_label + '_' + state] = true_a_probs[state]
                        scores_dict[args.group_B_label + '_' + state] = true_b_probs[state]

                if keep_max_scores_per_bin:
                    min_state = min(states, key=lambda s: scores_dict[s])
                    max_state = max(states, key=lambda s: scores_dict[s])
                    for s in states:
                        if s not in [min_state, max_state]:
                            scores_dict[s] = 0

            if not background_chunk:

                true_a_bins = set(a_bins if true_a_idx == -1 else a_window[true_a_idx * len(datasets_A):(true_a_idx + 1) * len(datasets_A)])
                true_b_bins = set(b_bins if true_b_idx == -1 else b_window[true_b_idx * len(datasets_B):(true_b_idx + 1) * len(datasets_B)])

                if use_posteriors:
                    state_indices = range(len(states))
                    true_a_bins = set(states[max(state_indices, key=lambda i: p[i])] for p in true_a_bins)
                    true_b_bins = set(states[max(state_indices, key=lambda i: p[i])] for p in true_b_bins)

                    a_bins = [states[max(state_indices, key=lambda i: p[i])] for p in a_bins]
                    b_bins = [states[max(state_indices, key=lambda i: p[i])] for p in b_bins]

                closest_state_A = min(true_a_bins, key=lambda s: KL(true_a_probs, consistent_state_probs[GROUP_A][s]))
                closest_state_B = min(true_b_bins, key=lambda s: KL(true_b_probs, consistent_state_probs[GROUP_B][s]))

                closest_transition_cache[window_idx] = (closest_state_A, closest_state_B, a_bins, b_bins)

            scores_cache[window_idx] = scores_dict

        else:

            scores_dict = scores_cache[window_idx]

        # compute the score between A and B
        scores[bin_idx] = scores_dict

        if not background_chunk:
            state_transitions[bin_idx] = closest_transition_cache[window_idx]

    if verbose:
        echo(chrom_length, len(scores_cache), 100 * (cache_miss / float(chrom_length)), cache_miss)

    return expand_non_overlapping_windows(chrom, scores, state_transitions, chrom_segmentations_windows, background_chunk)



def store_wig(chrom,
              signal,
              out_signal_f,
              span):

    out_signal_f.write('fixedStep\tchrom=%s\tstart=0\tstep=%d\tspan=%d\n' % (chrom, span, span))

    for bin_signal in signal:

        out_signal_f.write('%.5lf\n' % bin_signal)


def store_dict_wig(chrom,
                   signal,
                   out_signal_files,
                   span):

    for signal_type in signal:
        # echo(signal_type)
        store_wig(chrom, signal[signal_type], out_signal_files[signal_type], span)


import bisect

def get_FDR(signal, background_model, cache):

    if signal not in cache:
        idx = bisect.bisect_left(background_model, (signal, 0))

        if idx == len(background_model):
            cache[signal] = background_model[-1][1]
        else:
            bm_signal, bm_fdr = background_model[idx]

            if signal == bm_signal:
                cache[signal] = bm_fdr
            else:
                # interpolate the FDR
                prev_signal, prev_fdr = background_model[idx - 1]
                cache[signal] = prev_fdr
                # get_FDR_cache[key] = ((bm_fdr - prev_fdr) / (bm_signal - prev_signal)) * signal + (prev_fdr * bm_signal - bm_fdr * prev_signal) / (bm_signal - prev_signal)

    return cache[signal]


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
                                                         states,

                                                         metric_A,
                                                         metric_B,

                                                         n_perm=100,
                                                         to_smooth=False,
                                                         max_min_window=0,
                                                         compute_per_state_scores=False,
                                                         n_threads=1,
                                                         posteriors_dir=None,
                                                         use_mean_distance_matrix=False,
                                                         keep_max_scores_per_bin=False
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
        n_comb = n_combinations(len(all_seg_keys), n_A_segs) / 2
    else:
        n_comb = n_combinations(len(all_seg_keys), n_A_segs)
    echo('Total combos:', n_comb)

    def shuffled_samples(n_perm):

        seen = set()

        if n_comb > n_perm:

            # yield A_segs_keys, B_segs_keys
            # seen = set([(A_segs_keys, B_segs_keys),
            #             (B_segs_keys, A_segs_keys)])
            # n_perm -= 1

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
        import functools
        # _map = functools.partial(pool.imap_unordered, chunksize=(n_perm / n_threads))
        _map = functools.partial(pool.imap_unordered, chunksize=max(1, (n_perm / n_threads)))

    i = 0
    # all_metrics = dict((k, v) for m in [metric_A, metric_B] for k, v in m.iteritems())

    for combo_background_scores in _map(worker,
                                        ((process_shuffled_combo,
                                          combo_no,

                                          # dict((d, all_metrics[d]) for d in shuffled_A_datasets),
                                          shuffled_A_datasets,
                                          dict((dataset, all_segs[dataset]) for dataset in shuffled_A_datasets),

                                          # dict((d, all_metrics[d]) for d in shuffled_B_datasets),
                                          shuffled_B_datasets,
                                          dict((dataset, all_segs[dataset]) for dataset in shuffled_B_datasets),

                                          BIN_SIZE,
                                          states,
                                          compute_per_state_scores,
                                          max_min_window,
                                          to_smooth,
                                          posteriors_dir,
                                          use_mean_distance_matrix,
                                          keep_max_scores_per_bin
                                          )
                                         for combo_no, (shuffled_A_datasets, shuffled_B_datasets) in enumerate(shuffled_samples(n_perm)))):

        for score_type in combo_background_scores:
            for s in combo_background_scores[score_type]:
                if s not in background_scores[score_type]:
                    background_scores[score_type][s] = 0
                background_scores[score_type][s] += combo_background_scores[score_type][s]
        i += 1
        echo('Main thread update:', i)
        gc.collect()

    if n_threads > 1:
        pool.close()

    for score_type in background_scores:
        total_regions = sum(background_scores[score_type].itervalues())
        total_so_far = 0

        for s in sorted(background_scores[score_type], reverse=True):
            total_so_far += background_scores[score_type][s]
            background_scores[score_type][s] = total_so_far / float(total_regions)

    clear_metric_cache()

    return dict((score_type, sorted(background_scores[score_type].items())) for score_type in background_scores)


def worker(args):
    func = None

    try:
        func = args[0]
        return func(*args[1:])

    except Exception, e:
        print 'Caught exception in output worker thread (pid: %d):' % os.getpid()
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

                            shuffled_B_datasets,
                            shuffled_B_segmentations,

                            BIN_SIZE,
                            states,
                            compute_per_state_scores,
                            max_min_window,
                            to_smooth,
                            posteriors_dir,
                            use_mean_distance_matrix,
                            keep_max_scores_per_bin):

    echo('Combo no:', combo_no)

    echo('Shuffled A:', shuffled_A_datasets)
    echo('Shuffled B:', shuffled_B_datasets)
    print

    shuffled_metric_A = learn_metric_from_all_replicates(shuffled_A_segmentations, states, BIN_SIZE, posteriors_dir)
    # learn_metric(shuffled_A_segmentations, states, BIN_SIZE)
    # print_metric(shuffled_metric_A)
    #
    shuffled_metric_B = learn_metric_from_all_replicates(shuffled_B_segmentations, states, BIN_SIZE, posteriors_dir)

    # print_metric(shuffled_metric_B)
    # print '*' * 50
    # print_average_metric(states, shuffled_metric_A, shuffled_metric_B)
    # print '*' * 50
    if use_mean_distance_matrix:
        shuffled_metric_A, shuffled_metric_B = compute_average_metric(shuffled_metric_A, shuffled_metric_B, states)

    # iterate over the bins in both segmentation groups
    # for chrom in random.sample(chromosomes, 2):

    background_scores = dict((k, {}) for k in [OVERALL_SCORE] + (states if compute_per_state_scores else []))
    random_chunks = generate_random_chunks(shuffled_A_segmentations, chromosomes, BIN_SIZE, N_RANDOM_CHUNKS=10)

    for chrom in sorted(random_chunks):
        # echo(combo_no, chrom)
        if posteriors_dir is None:
            chrom_segmentations_A = dict((d, chrom_segmentation_to_list(shuffled_A_segmentations[d][chrom], BIN_SIZE)) for d in shuffled_A_datasets)
            chrom_segmentations_B = dict((d, chrom_segmentation_to_list(shuffled_B_segmentations[d][chrom], BIN_SIZE)) for d in shuffled_B_datasets)
        else:
            chrom_segmentations_A = dict((d, read_posteriors_as_list(get_chrom_posteriors_fname(posteriors_dir, d, chrom))) for d in shuffled_A_datasets)
            chrom_segmentations_B = dict((d, read_posteriors_as_list(get_chrom_posteriors_fname(posteriors_dir, d, chrom))) for d in shuffled_B_datasets)

        for chunk_start, chunk_end in random_chunks[chrom]:
            # echo(combo_no, chunk_start, chunk_end)

            chunk_scores = get_overall_and_per_state_diff_score(chrom,
                                                                slice_segmentations(chrom_segmentations_A, chunk_start, chunk_end),
                                                                shuffled_metric_A,
                                                                slice_segmentations(chrom_segmentations_B, chunk_start, chunk_end),
                                                                shuffled_metric_B,
                                                                states,
                                                                None,
                                                                compute_per_state_scores,
                                                                max_min_window,
                                                                background_chunk=True,
                                                                use_posteriors=(posteriors_dir is not None),
                                                                keep_max_scores_per_bin=keep_max_scores_per_bin)

            if to_smooth:
                smooth_dict(chunk_scores) #chunk_scores = dict((score_type, smooth(chunk_scores[score_type])) for score_type in chunk_scores)

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
                              states,
                              metric_A,
                              metric_B,
                              to_smooth=False,
                              max_min_window=0,
                              compute_per_state_scores=False,
                              posteriors_dir=None,
                              keep_max_scores_per_bin=False):

    echo('Computing foreground scores')
    foreground_scores = dict((k, {}) for k in [OVERALL_SCORE] + (states if compute_per_state_scores else []))

    datasets_A = sorted(group_A_segmentations)
    datasets_B = sorted(group_B_segmentations)

    random_chunks = generate_random_chunks(group_A_segmentations, chromosomes, BIN_SIZE, N_RANDOM_CHUNKS=100)

    for chrom in random_chunks:

        if posteriors_dir is None:
            chrom_segmentations_A = dict((d, chrom_segmentation_to_list(group_A_segmentations[d][chrom], BIN_SIZE)) for d in datasets_A)
            chrom_segmentations_B = dict((d, chrom_segmentation_to_list(group_B_segmentations[d][chrom], BIN_SIZE)) for d in datasets_B)

        else:
            chrom_segmentations_A = dict((d, read_posteriors_as_list(get_chrom_posteriors_fname(posteriors_dir, d, chrom))) for d in datasets_A)
            chrom_segmentations_B = dict((d, read_posteriors_as_list(get_chrom_posteriors_fname(posteriors_dir, d, chrom))) for d in datasets_B)

        for chunk_no, (chunk_start, chunk_end) in enumerate(random_chunks[chrom]):

            chunk_scores = get_overall_and_per_state_diff_score(chrom + '_' + str(chunk_no),
                                                                slice_segmentations(chrom_segmentations_A, chunk_start, chunk_end),
                                                                metric_A,
                                                                slice_segmentations(chrom_segmentations_B, chunk_start, chunk_end),
                                                                metric_B,
                                                                states,
                                                                None,
                                                                compute_per_state_scores,
                                                                max_min_window,
                                                                background_chunk=True,
                                                                use_posteriors=(posteriors_dir is not None),
                                                                keep_max_scores_per_bin=keep_max_scores_per_bin)

            if to_smooth:
                smooth_dict(chunk_scores) #chunk_scores = dict((score_type, smooth(chunk_scores[score_type])) for score_type in chunk_scores)

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
                             group_A_segmentations,
                             metric_A,
                             group_B_segmentations,
                             metric_B,
                             states,
                             posteriors_dir=None,
                             use_mean_distance_matrix=False,
                             keep_max_scores_per_bin=False):

    # pick the random chunks
    BIN_SIZE = args.bin_size
    MAX_N_SHUFFLED_SAMPLES = 100

    if args.chrom_hmm_binarized:
        from permute_marks_in_binarized_ChromHMM import compute_background_scores_by_shuffling_marks
        seg = group_A_segmentations.values()[0]
        chrom_lengths = dict((chrom, seg[chrom][-1][1] / BIN_SIZE) for chrom in seg)

        # print 'chrom lengths:'
        # print chrom_lengths

        echo('Learning significance threshold by permuting histone marks at q-value:', fdr_threshold)

        background_scores = compute_background_scores_by_shuffling_marks(args.chrom_hmm_binarized,
                                                                         group_A_segmentations.keys() + group_B_segmentations.keys(),
                                                                         args.chrom_hmm_model_path,
                                                                         chrom_lengths,

                                                                         n_group_A=len(group_A_segmentations),
                                                                         BIN_SIZE=BIN_SIZE,
                                                                         states=states,
                                                                         n_perms=MAX_N_SHUFFLED_SAMPLES,
                                                                         compute_per_state_scores=args.per_state_scores,
                                                                         max_min_window=args.max_min_window,
                                                                         to_smooth=args.smooth,
                                                                         use_posteriors=(posteriors_dir is not None),
                                                                         n_threads=args.n_threads,
                                                                         use_mean_distance_matrix=use_mean_distance_matrix,
                                                                         keep_max_scores_per_bin=keep_max_scores_per_bin)

    else:

        echo('Learning significance threshold by permuting segmentations at q-value:', fdr_threshold)

        background_scores = compute_background_scores_by_shuffling_segmentations(group_A_segmentations,
                                                                                 group_B_segmentations,
                                                                                 states,

                                                                                 metric_A,
                                                                                 metric_B,

                                                                                 n_perm=MAX_N_SHUFFLED_SAMPLES,
                                                                                 compute_per_state_scores=args.per_state_scores,
                                                                                 to_smooth=args.smooth,
                                                                                 max_min_window=args.max_min_window,
                                                                                 n_threads=args.n_threads,
                                                                                 posteriors_dir=posteriors_dir,
                                                                                 use_mean_distance_matrix=use_mean_distance_matrix,
                                                                                 keep_max_scores_per_bin=keep_max_scores_per_bin)

    foreground_scores = compute_foreground_scores(group_A_segmentations,
                                                  group_B_segmentations,
                                                  states,
                                                  metric_A,
                                                  metric_B,
                                                  to_smooth=args.smooth,
                                                  max_min_window=args.max_min_window,
                                                  compute_per_state_scores=args.per_state_scores,
                                                  posteriors_dir=posteriors_dir,
                                                  keep_max_scores_per_bin=keep_max_scores_per_bin)

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
            # background_threshold[score_type] = all_scores[-1][0] + 1

            # background_threshold[score_type] = best_score_cutoff

    return dict((score_type, sorted(background_model[score_type].iteritems())) for score_type in background_model), background_threshold

SIGNAL = 'signal'
STATE_TRANSITIONS = 'state_transitions'

def process_chromosome(out_prefix,

                       chrom,
                       group_A_segmentations,
                       metric_A,
                       group_B_segmentations,
                       metric_B,
                       states,
                       consistent_state_probs,
                       compute_per_state_scores,
                       max_min_window,
                       background_chunk=False,
                       use_posteriors=False,
                       report_per_group_state_probabilities=False,
                       keep_max_scores_per_bin=False):

    echo('Processing:', chrom)
    chrom, signal, state_transitions = get_overall_and_per_state_diff_score(chrom,
                                                                            group_A_segmentations,
                                                                            metric_A,
                                                                            group_B_segmentations,
                                                                            metric_B,
                                                                            states,
                                                                            consistent_state_probs,
                                                                            compute_per_state_scores,
                                                                            max_min_window,
                                                                            background_chunk,
                                                                            use_posteriors,
                                                                            report_per_group_state_probabilities,
                                                                            keep_max_scores_per_bin)

    if args.smooth:
        echo('Smoothing:', chrom)
        smooth_dict(signal)

    echo('Writing down pickle files:', chrom)

    st_fname = out_prefix + '.' + chrom + '.' + STATE_TRANSITIONS + '.pickle'
    with open(st_fname, 'w') as out_f:
        pickle.dump(state_transitions, out_f, pickle.HIGHEST_PROTOCOL)

    pickle_fnames = {}
    for score_type in signal:
        pickle_fname = out_prefix + '.' + chrom + '.' + score_type + '.pickle'
        pickle_fnames[score_type] = pickle_fname

        with open(pickle_fname, 'w') as out_f:
            pickle.dump(signal[score_type], out_f, pickle.HIGHEST_PROTOCOL)

    echo('Done:', chrom)

    return pickle_fnames, st_fname, chrom


def merge_output_files(out_prefix,
                       score_type,
                       score_fnames,
                       st_fnames,
                       BIN_SIZE,
                       score_background_model,
                       threshold,
                       skip_bed_files=False):

    echo('Merging:', score_type)
    wig_fname = out_prefix + '.' + re.sub(r'\W+', '_', score_type) + '.wig.gz'
    bed_fname = out_prefix + '.' + re.sub(r'\W+', '_', score_type) + '.summary.bed.gz'

    for fname in [wig_fname, bed_fname]:
        if os.path.isfile(fname):
            os.unlink(fname)

    echo('Reading pickle files:', score_type)

    signal = {}
    state_transitions = {}

    for chrom in score_fnames:
        with open(score_fnames[chrom]) as in_f:
            signal[chrom] = pickle.load(in_f)
        os.unlink(score_fnames[chrom])

        with open(st_fnames[chrom]) as in_f:
            state_transitions[chrom] = pickle.load(in_f)

    echo('Writing down wig file:', wig_fname)
    with open_file(wig_fname, 'w') as out_f:
        title = os.path.split(out_prefix)[1].replace('.wig', '').replace('.gz', '') + ' ' + score_type
        out_f.write('track type=wiggle_0 name="%s" description="%s"\n' % (title, title))
        for chrom in sorted(signal):
            store_wig(chrom, signal[chrom], out_f, BIN_SIZE)

    if skip_bed_files:
        return True

    echo('Sorting scores')

    bins_sorted_by_score = sorted([(chrom, bin_no)
                                    for chrom in signal
                                        for bin_no in xrange(len(signal[chrom]))],
                           key=lambda (chrom, bin_no): (-signal[chrom][bin_no], chrom, bin_no))

    # sorted_scores = sorted([(score, chrom, bin_no, state_transition)
    #                                 for chrom in signal
    #                                     for bin_no, (score, state_transition) in enumerate(izip(signal[chrom],
    #                                                                                             state_transitions[chrom]))],
    #                        key=lambda (score, chrom, bin_no, state_transition): (-score, chrom, bin_no))
    #
    echo('Writing sorted summary files for', score_type, ':', bed_fname)

    score_cache = {}

    with open_file(bed_fname, 'w') as out_f:
        # for score, chrom, bin_no, (closest_state_A, closest_state_B, a_states, b_states) in sorted_scores:
        for chrom, bin_no in bins_sorted_by_score:
            score = signal[chrom][bin_no]
            (closest_state_A, closest_state_B, a_states, b_states) = state_transitions[chrom][bin_no]

            out_f.write('%s\t%d\t%d\t.\t%lf\t.\t%s\t%s\t%lf\t%s\n' %
                        (chrom,
                         bin_no * BIN_SIZE,
                         (bin_no + 1) * BIN_SIZE,
                         score,

                         closest_state_A + '-' + closest_state_B,
                         ','.join(a_states) + '/' + ','.join(b_states),
                         get_FDR(abs(score), score_background_model, score_cache),
                         '1' if threshold is not None and abs(score) >= threshold else '0'))

    gc.collect()

    echo('Done:', score_type)
    return True

# def merge_output_files(out_prefix,
#                        score_type,
#                        score_fnames,
#                        st_fnames,
#                        BIN_SIZE,
#                        score_background_model,
#                        threshold,
#                        skip_bed_files=False):
#
#     echo('Merging:', score_type)
#     wig_fname = out_prefix + '.' + re.sub(r'\W+', '_', score_type) + '.wig.gz'
#     bed_fname = out_prefix + '.' + re.sub(r'\W+', '_', score_type) + '.summary.bed.gz'
#
#     for fname in [wig_fname, bed_fname]:
#         if os.path.isfile(fname):
#             os.unlink(fname)
#
#     echo('Reading pickle files:', score_type)
#
#     signal = {}
#     state_transitions = {}
#
#     for chrom in score_fnames:
#         with open(score_fnames[chrom]) as in_f:
#             signal[chrom] = pickle.load(in_f)
#         os.unlink(score_fnames[chrom])
#
#         with open(st_fnames[chrom]) as in_f:
#             state_transitions[chrom] = pickle.load(in_f)
#
#     echo('Writing down wig file:', wig_fname)
#     with open_file(wig_fname, 'w') as out_f:
#         title = os.path.split(out_prefix)[1].replace('.wig', '').replace('.gz', '') + ' ' + score_type
#         out_f.write('track type=wiggle_0 name="%s" description="%s"\n' % (title, title))
#         for chrom in sorted(signal):
#             store_wig(chrom, signal[chrom], out_f, BIN_SIZE)
#
#     if skip_bed_files:
#         return True
#
#     echo('Sorting scores')
#     sorted_scores_indexes = dict((c, 0) for c in signal)
#     sorted_scores = {}
#     for chrom in signal.keys():
#         sorted_scores[chrom] = sorted([(score, bin_no, state_transition)
#                                             for bin_no, (score, state_transition) in enumerate(izip(signal[chrom],
#                                                                                                     state_transitions[chrom]))],
#                                       key=lambda (score, bin_no, state_transition): (-score, bin_no))
#
#         del signal[chrom]
#         del state_transitions[chrom]
#     gc.collect()
#     # sorted_scores = sorted([(score, chrom, bin_no, state_transition)
#     #                                 for chrom in signal
#     #                                     for bin_no, (score, state_transition) in enumerate(izip(signal[chrom],
#     #                                                                                             state_transitions[chrom]))],
#     #                        key=lambda (score, chrom, bin_no, state_transition): (-score, chrom, bin_no))
#
#     echo('Writing sorted summary files for', score_type, ':', bed_fname)
#
#     score_cache = {}
#     def get_score(chrom, sorted_scores, sorted_scores_indexes):
#         idx = sorted_scores_indexes[chrom]
#
#         if idx == len(sorted_scores[chrom]):
#             return None
#
#         sorted_scores_indexes[chrom] += 1
#         return sorted_scores[chrom][idx], chrom
#
#     chromosomes = sorted(sorted_scores)
#
#     c_scores = dict((chrom, get_score(chrom, sorted_scores, sorted_scores_indexes)) for chrom in chromosomes)
#
#     with open_file(bed_fname, 'w') as out_f:
#         while len(c_scores) > 0:
#             (score, bin_no, (closest_state_A, closest_state_B, a_states, b_states)), chrom = max(c_scores.itervalues())
#
#         # for score, chrom, bin_no, (closest_state_A, closest_state_B, a_states, b_states) in sorted_scores:
#             out_f.write('%s\t%d\t%d\t.\t%lf\t.\t%s\t%s\t%lf\t%s\n' %
#                         (chrom,
#                          bin_no * BIN_SIZE,
#                          (bin_no + 1) * BIN_SIZE,
#                          score,
#
#                          closest_state_A + '-' + closest_state_B,
#                          ','.join(a_states) + '/' + ','.join(b_states),
#                          get_FDR(abs(score), score_background_model, score_cache),
#                          '1' if threshold is not None and abs(score) >= threshold else '0'))
#
#             new_score = get_score(chrom, sorted_scores, sorted_scores_indexes)
#             if new_score is None:
#                 del c_scores[chrom]
#             else:
#                 c_scores[chrom] = new_score
#
#     gc.collect()
#
#     echo('Done:', score_type)
#     return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-a', nargs='+', help='segmentation files of group A')
    parser.add_argument('-b', nargs='+', help='segmentation files of group B')

    parser.add_argument('-A', dest='A', help='text file with segmentation files of group A, one per line')
    parser.add_argument('-B', dest='B', help='text file with segmentation files of group B, one per line')

    parser.add_argument('--bin-size', type=int, help='Bin size, default: %(default)s', default=200)
    parser.add_argument('-f', '--fdr-threshold', type=float, help='FDR threshold: %(default)s', default=0.05)

    parser.add_argument('--group-A-label', help='label for group A', default=GROUP_A)
    parser.add_argument('--group-B-label', help='label for group B', default=GROUP_B)

    parser.add_argument('-o', '--output', help='output prefix')
    # parser.add_argument('--background-scores', help='bacground_scores.pickle')

    parser.add_argument('--chrom-hmm-binarized', help='path to ChromHMM binarized files to compute background scores')
    parser.add_argument('--chrom-hmm-model-path', help='path to the ChromHMM model to use for shuffled marks')

    parser.add_argument('--posteriors-dir', dest='posteriors_dir', help='posteriors directory from ChromHMM '
                                                                        '(scores will be computed by '
                                                                        'taking into account the posteriors')

    parser.add_argument('--background-model', dest='background_model', help='Pickle file with previously '
                                                                            'computed background model')

    parser.add_argument('-s', '--smooth', action='store_true', help='smooth signal')
    parser.add_argument('-p', '--n-threads', type=int, help='number of threads to use', default=1)

    parser.add_argument('--max-min-window', type=int, help='pick the maximum distance between the two most similar bins within this window', default=0)

    parser.add_argument('--filter-chromosomes', nargs='+', dest='filter_chroms',
                        help='Apply method only on these chromosomes. Default: use all chromosomes',
                        default=None)

    parser.add_argument('--keep-max-scores-per-bin',
                        dest='keep_max_scores_per_bin',
                        action='store_true',
                        help='Keep only the maximum score in each direction per bin. Default: %(default)s',
                        default=False)

    parser.add_argument('--per-state-scores', action='store_true', help='compute per state scores', default=False)

    parser.add_argument('--report-per-group-state-probabilities',
                        action='store_true',
                        help='output the probabilities for each state in each group per state scores',
                        default=False)
    parser.add_argument('--skip-bed-files', action='store_true', help='output only wiggle files', default=False)
    parser.add_argument('--use-mean-distance-matrix', action='store_true',
                        help='Use one distance matrix for all samples computed as the average distance between each '
                             'pair of states',
                        default=False)

    # parser.add_argument('--use-closest-rep', action='store_true', help='use only closest replicate to learn the metric')
    # parser.add_argument('--all-for-null', action='store_true', default=False, help='include the original combo in the background model')
    # parser.add_argument('--store-scores', action='store_true', default=False, help='store foreground and background scores in a pickle file')
    # parser.add_argument('--output-changes', action='store_true', default=False, help='output significant pairwise changes')

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)

    open_log(args.output + '.log')
    BIN_SIZE = args.bin_size
    echo('cmd:', ' '.join(sys.argv))
    for arg in sorted(vars(args)):
        echo(arg, '=', getattr(args, arg))

    keep_max_scores_per_bin = args.keep_max_scores_per_bin
    posteriors_dir = args.posteriors_dir
    use_mean_distance_matrix = args.use_mean_distance_matrix

    max_min_window = args.max_min_window

    compute_per_state_scores = args.per_state_scores

    if args.A:
        seg_A_fnames = [os.path.join(os.path.split(args.A)[0], l.strip()) for l in open_file(args.A)]
    else:
        seg_A_fnames = args.a

    if args.B:
        seg_B_fnames = [os.path.join(os.path.split(args.B)[0], l.strip()) for l in open_file(args.B)]
    else:
        seg_B_fnames = args.b

    group_A_segmentations, states_A = read_segmentations(seg_A_fnames)
    group_B_segmentations, states_B = read_segmentations(seg_B_fnames)

    states = sorted(set(states_A) | set(states_B), key=state_key)

    if args.filter_chroms is not None:
        chromosomes = args.filter_chroms
        echo('Using chromosomes:', chromosomes)
    else:
        chromosomes = sorted(reduce(operator.and_, [set(s.keys())
                                                 for segs in [group_A_segmentations, group_B_segmentations]
                                                    for s in segs.values()]),
                                    key=lambda c: int(c.replace('chr','')) if re.search(r'^chr(\d+)$', c) else 100)
    # print chromosomes
    group_A_segmentations = filter_chroms(group_A_segmentations, chromosomes)
    group_B_segmentations = filter_chroms(group_B_segmentations, chromosomes)

    # print 'A'
    # metric_A = learn_metric_from_all_replicates_test(group_A_segmentations, states, BIN_SIZE, posteriors_dir)
    # print 'B'
    # metric_B = learn_metric_from_all_replicates_test(group_B_segmentations, states, BIN_SIZE, posteriors_dir)
    # exit(1)

    # posteriors_A = learn_posterior_from_all_replicates(group_A_segmentations, states, BIN_SIZE)
    # exit(1)
    metric_A = learn_metric_from_all_replicates(group_A_segmentations, states, BIN_SIZE, posteriors_dir)
    metric_B = learn_metric_from_all_replicates(group_B_segmentations, states, BIN_SIZE, posteriors_dir)

    if use_mean_distance_matrix:
        metric_A, metric_B = compute_average_metric(metric_A, metric_B, states)

    echo('Metric A')
    print_metric(metric_A)

    echo('Metric B')
    print_metric(metric_B)

    print_consistent_state_transition_scores(states, metric_A, metric_B)

    _average_metric, _ = compute_average_metric(metric_A, metric_B, states)
    echo('Average metric')
    print_metric(_average_metric)

    fdr_threshold = args.fdr_threshold

    out_fname = args.output

    score_types = [OVERALL_SCORE] + (states if compute_per_state_scores else [])
    if args.fdr_threshold < 1:
        if args.background_model is not None:
            background_model, background_threshold = pickle.load(open(args.background_model))
        else:
            background_model, background_threshold = compute_background_model(args,
                                                                              fdr_threshold,
                                                                              group_A_segmentations,
                                                                              metric_A,
                                                                              group_B_segmentations,
                                                                              metric_B,
                                                                              states,
                                                                              posteriors_dir,
                                                                              use_mean_distance_matrix,
                                                                              keep_max_scores_per_bin)

            with open(out_fname + '.background_model.pickle', 'w') as out_f:
                echo('Saving background model in:', out_fname + '.background_model.pickle')
                pickle.dump((background_model, background_threshold), out_f, pickle.HIGHEST_PROTOCOL)

        gc.collect()
    else:
        echo('Skipping null model computations')
        background_model = dict((s, [(0, 0)]) for s in score_types)
        background_threshold = dict((s, 0) for s in score_types)

    real_scores = {}
    real_state_transitions = {}
    echo('Computing real scores')

    consistent_state_probs = generate_consistent_state_probabilities(states, metric_A, metric_B)
    # chromosomes = chromosomes # ['chr1', 'chr2', 'chrM', 'chrY']
    if args.n_threads > 1:
        pool = Pool(args.n_threads)
        _map = pool.imap_unordered

    score_fnames = {}
    st_fnames = {}

    for pickle_fnames, st_fname, chrom in _map(worker,
                                                 ((process_chromosome,

                                                   out_fname,
                                                   _chrom,
                                                   dict((d, chrom_segmentation_to_list(group_A_segmentations[d][_chrom], BIN_SIZE)
                                                                 if posteriors_dir is None else
                                                            read_posteriors_as_list(get_chrom_posteriors_fname(posteriors_dir, d, _chrom))
                                                         )
                                                        for d in group_A_segmentations),
                                                   metric_A,
                                                   dict((d, chrom_segmentation_to_list(group_B_segmentations[d][_chrom], BIN_SIZE)
                                                                 if posteriors_dir is None else
                                                            read_posteriors_as_list(get_chrom_posteriors_fname(posteriors_dir, d, _chrom))
                                                        )
                                                        for d in group_B_segmentations),
                                                   metric_B,
                                                   states,
                                                   consistent_state_probs,
                                                   compute_per_state_scores,
                                                   max_min_window,
                                                   False,
                                                   (posteriors_dir is not None),
                                                   args.report_per_group_state_probabilities,
                                                   keep_max_scores_per_bin) for _chrom in chromosomes)):

        st_fnames[chrom] = st_fname
        for score_type in pickle_fnames:
            if score_type not in score_fnames:
                score_fnames[score_type] = {}

            score_fnames[score_type][chrom] = pickle_fnames[score_type]


    all(_map(worker, ((merge_output_files,
                       out_fname,
                       score_type,
                       score_fnames[score_type],
                       st_fnames,
                       BIN_SIZE,
                       background_model.get(score_type, [(0, 0)]),
                       background_threshold.get(score_type, None),
                       args.skip_bed_files
                       ) for score_type in score_fnames)))

    # delete pickled state transition files
    for fname in st_fnames.values():
        os.unlink(fname)

    if args.n_threads > 1:
        pool.close()

    echo('Done')
    close_log()
