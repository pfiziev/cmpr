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
                                         group_segmentations,
                                         metric,
                                         states,
                                         consistent_state_probs,
                                         compute_per_state_scores,
                                         max_min_window,
                                         background_chunk=False,
                                         use_posteriors=False,
                                         report_per_group_state_probabilities=False,
                                         constitutive_only=False,
                                         compute_differential_from_closest_group=False):

    echo('Chromosome:', chrom)

    groups = sorted(group_segmentations)

    # echo('Groups:', groups)
    datasets = dict((g, sorted(group_segmentations[g])) for g in groups)

    chrom_segmentations = dict((g, [group_segmentations[g][d] for d in datasets[g]]) for g in groups)

    chrom_length = len(chrom_segmentations.values()[0][0])
    echo('Chrom length:', chrom_length)

    score_types = [OVERALL_SCORE] + (states if compute_per_state_scores else [])
    if constitutive_only:
        scores = {}
    else:
        scores = dict((g, dict((k, [0.0] * chrom_length) for k in score_types)) for g in groups)

    scores[CONSTITUTIVE] = dict((k, [0.0] * chrom_length) for k in states)

    if constitutive_only:
        closest_state = {}
    else:
        closest_state = dict((g, [None] * chrom_length) for g in groups)

    closest_state[CONSTITUTIVE] = [None] * chrom_length

    # return scores, state_transitions
    if use_posteriors:
        _get_prob_vector = get_prob_vector_based_on_posteriors
    else:
        _get_prob_vector = get_prob_vector

    all_summed_probs = dict((s, 0) for s in states)
    rest_summed_probs = dict((s, 0) for s in states)

    n_groups = len(groups)

    group_to_group_distance = [[0] * n_groups for _ in xrange(n_groups)]

    for bin_idx in xrange(chrom_length):
        if bin_idx % 20000 == 0:
            echo('Processed bins:', bin_idx, 'on', chrom)

        all_current_bins = [tuple(seg[bin_idx] for seg in chrom_segmentations[g]) for g in groups]

        prob_vectors = [_get_prob_vector(bins, datasets[g], metric[g], states) for bins, g in izip(all_current_bins, groups)]

        for s in states:
            all_summed_probs[s] = sum(p[s] for p in prob_vectors)

        if compute_per_state_scores:
            for s in states:
                scores[CONSTITUTIVE][s][bin_idx] = all_summed_probs[s]

            # closest_state[CONSTITUTIVE][bin_idx] = (min([(b, i) for i in xrange(n_groups)
            #                                             for b in all_current_bins[i]],
            #                                                  key=lambda (s, i): KL(prob_vectors[i], consistent_state_probs[groups[i]][s]))[0],
            #
            #                                             [max(bins, key=lambda s: prob_vectors[i][s])
            #                                                 for i, bins in enumerate(all_current_bins)])

            # closest_state[CONSTITUTIVE][bin_idx] = min([(b, i) for i in xrange(n_groups)
            #                                                 for b in all_current_bins[i]],
            #                                                      key=lambda (s, i): KL(prob_vectors[i], consistent_state_probs[groups[i]][s]))[0]

            closest_state[CONSTITUTIVE][bin_idx] = max(states, key=lambda s: scores[CONSTITUTIVE][s][bin_idx])


        if constitutive_only:
            continue

        if compute_differential_from_closest_group:
            for g_idx_1 in xrange(n_groups):
                for g_idx_2 in xrange(g_idx_1 + 1, n_groups):
                    sym_kl_div = symmetric_KL_divergence(prob_vectors[g_idx_1], prob_vectors[g_idx_2])

                    group_to_group_distance[g_idx_1][g_idx_2] = sym_kl_div
                    group_to_group_distance[g_idx_2][g_idx_1] = sym_kl_div

        for g_idx_1, g, g_prob, g_bins in izip(xrange(n_groups), groups, prob_vectors, all_current_bins):

            if report_per_group_state_probabilities:
                for s in states:
                    scores[g][s][bin_idx] = g_prob[s]

            else:
                if compute_differential_from_closest_group:
                    closest_group_idx = min([g_idx_2 for g_idx_2 in xrange(n_groups) if g_idx_2 != g_idx_1],
                                             key=lambda g_idx_2: group_to_group_distance[g_idx_1][g_idx_2])

                    scores[g][OVERALL_SCORE][bin_idx] = symmetric_KL_divergence(g_prob, prob_vectors[closest_group_idx])

                    if compute_per_state_scores:
                        for s in states:
                            scores[g][s][bin_idx] = g_prob[s] - prob_vectors[closest_group_idx][s]

                else:
                    for s in states:
                        rest_summed_probs[s] = (all_summed_probs[s] - g_prob[s]) / (n_groups - 1)

                    scores[g][OVERALL_SCORE][bin_idx] = symmetric_KL_divergence(g_prob, rest_summed_probs)
                    if compute_per_state_scores:
                        for s in states:
                            scores[g][s][bin_idx] = g_prob[s] - rest_summed_probs[s]

            # closest_state[g][bin_idx] = (min(g_bins, key=lambda s: KL(g_prob, consistent_state_probs[g][s])), g_bins)
            closest_state[g][bin_idx] = min(g_bins, key=lambda s: KL(g_prob, consistent_state_probs[g][s]))
    # echo('Done')
    # exit(1)
    return chrom, scores, closest_state



import ctypes
clib = ctypes.CDLL(os.path.join(os.path.split(__file__)[0],
                                    'csdelta.so'))
C_DOUBLE = ctypes.c_double

def get_overall_and_per_state_diff_score_C(chrom,
                                           C_chrom_length,
                                           group_segmentations_compressed,
                                           metric,
                                           states,
                                           consistent_state_probs,
                                           compute_per_state_scores,
                                           max_min_window,
                                           background_chunk=False,
                                           use_posteriors=False,
                                           report_per_group_state_probabilities=False,
                                           constitutive_only=False,
                                           compute_differential_from_closest_group=False):

    echo('Chromosome:', chrom)


    groups = sorted(group_segmentations_compressed)
    # echo('Groups:', groups)

    datasets = dict((g, sorted(group_segmentations_compressed[g])) for g in groups)
    # for g_idx, g in enumerate(groups):
    #     echo(g_idx, 'Group:', g, 'Datasets:', datasets[g])

    n_datasets = sum(len(datasets[g]) for g in groups)

    C_n_groups = len(groups)
    C_n_states = len(states)

    chrom_segmentations = dict((g, [chrom_segmentation_to_list(decompress_chrom(group_segmentations_compressed[g][d]), BIN_SIZE)
                                    for d in datasets[g]]) for g in groups)

    state_to_int = dict((s, states.index(s)) for s in states)

    C_group_segmentations = (ctypes.c_int * (n_datasets * C_chrom_length))(*[
                                        state_to_int[dataset_segmentation[bin_idx]]
                                            for group in groups
                                                for dataset_segmentation in chrom_segmentations[group]
                                                    for bin_idx in xrange(C_chrom_length)])

    # echo('Memory freed from decompressed chrom_segmentations')

    C_datasets_per_group = (ctypes.c_int * C_n_groups)(*[len(datasets[g]) for g in groups])

    group_offset = [0] * C_n_groups
    for g_idx, g in enumerate(groups):
        if g_idx == 0:
            continue

        group_offset[g_idx] = group_offset[g_idx - 1] + C_datasets_per_group[g_idx - 1]

    C_group_offset = (ctypes.c_int * C_n_groups)(*group_offset)

    C_metric = (C_DOUBLE * (n_datasets * C_n_states * C_n_states))(*[
        metric[g][d][s1][s2]
        for g in groups
            for d in datasets[g]
                for s1 in states
                    for s2 in states
    ])

    C_consistent_state_probs = (C_DOUBLE * (C_n_groups * C_n_states * C_n_states))(*[
        consistent_state_probs[g][s1][s2]
        for g in groups
                for s1 in states
                    for s2 in states
    ])

    C_overall_scores = (C_DOUBLE * (C_n_groups * C_chrom_length))()
    C_per_state_scores = (C_DOUBLE * (C_n_groups * C_chrom_length * C_n_states))()
    C_constitutive_scores = (C_DOUBLE * (C_n_states * C_chrom_length))()

    C_closest_state_per_group = (ctypes.c_int * (C_n_groups * C_chrom_length))()
    C_closest_constitutive_state = (ctypes.c_int * (C_chrom_length))()

    echo('Chrom length:', C_chrom_length)

    clib.get_overall_and_per_state_diff_score_for_multiple_groups(C_group_segmentations,
                                                                  C_n_groups,

                                                                  C_group_offset,
                                                                  C_datasets_per_group,

                                                                  C_chrom_length,

                                                                  C_metric,
                                                                  C_n_states,

                                                                  C_consistent_state_probs,

                                                                  int(compute_per_state_scores),
                                                                  int(compute_differential_from_closest_group),

                                                                  # results variables
                                                                  C_overall_scores,
                                                                  C_per_state_scores,

                                                                  C_constitutive_scores,

                                                                  C_closest_state_per_group,
                                                                  C_closest_constitutive_state
                                                                  )



    # echo('Groups:', groups)

    scores = dict((g, {}) for g in groups)
    closest_state = {}

    scores[CONSTITUTIVE] = {}
    for state_idx, state in enumerate(states):
        scores[CONSTITUTIVE][state] = C_constitutive_scores[state_idx * C_chrom_length: (state_idx + 1) * C_chrom_length]

    closest_state[CONSTITUTIVE] = [states[state_idx] for state_idx in C_closest_constitutive_state]

    for g_idx, g in enumerate(groups):
        closest_state[g] = [states[C_closest_state_per_group[bin_idx]]
                            for bin_idx in xrange(g_idx * C_chrom_length, (g_idx + 1) * C_chrom_length)]

        scores[g][OVERALL_SCORE] = C_overall_scores[g_idx * C_chrom_length: (g_idx + 1) * C_chrom_length]

        for state_idx, state in enumerate(states):
            offset = g_idx * C_n_states * C_chrom_length + state_idx * C_chrom_length
            scores[g][state] = C_per_state_scores[offset: offset + C_chrom_length]


    return chrom, scores, closest_state



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


def compute_background_scores_by_shuffling_segmentations(group_segmentations,
                                                         states,

                                                         metric,

                                                         n_perm=100,
                                                         to_smooth=False,
                                                         max_min_window=0,
                                                         n_threads=1,
                                                         use_mean_distance_matrix=False,
                                                         compute_differential_from_closest_group=False
                                                         ):

    # longest_chromosome = max(chromosomes, key=lambda c: group_A_segmentations.values()[0][c][-1][1])
    # print 'Longest chromosome:', longest_chromosome

    groups = sorted(group_segmentations)
    n_segs = [len(group_segmentations[g]) for g in groups]

    datasets = dict((g, sorted(group_segmentations[g])) for g in groups)

    # make the union of the two groups
    all_seg_keys = tuple(d for g in groups for d in datasets[g])
    all_segs = dict((d, group_segmentations[g][d]) for g in groups for d in datasets[g])

    n_comb = 1
    _n_datasets = sum(n_segs)

    for g_idx, g in enumerate(groups):
        n_comb *= n_combinations(n_segs[g_idx], _n_datasets)
        _n_datasets -= n_segs[g_idx]

    if all(n_segs[0] == ns for ns in n_segs):
        n_comb /= 2

    echo('Total combos:', n_comb)

    def shuffled_samples(n_perm):

        seen = set()

        if n_comb > n_perm:

            # yield A_segs_keys, B_segs_keys
            # seen = set([(A_segs_keys, B_segs_keys),
            #             (B_segs_keys, A_segs_keys)])
            # n_perm -= 1
            def generate_random_grouping():
                shuffled_datasets = list(all_seg_keys)
                random.shuffle(shuffled_datasets)
                g_offset = 0
                shuffled_groups = []

                for g_idx in xrange(len(groups)):
                    shuffled_groups.append(tuple(sorted(shuffled_datasets[g_offset: g_offset + n_segs[g_idx]])))
                    g_offset += n_segs[g_idx]

                return tuple(shuffled_groups)

            for _ in xrange(n_perm):

                shuffled_grouping = generate_random_grouping()

                while shuffled_grouping in seen or tuple(reversed(shuffled_grouping)) in seen:
                    shuffled_grouping = generate_random_grouping()

                seen.add(shuffled_grouping)
                seen.add(tuple(reversed(shuffled_grouping)))

                yield shuffled_grouping

        else:
            # generate all shuffled samples

            for shuffled_A_segs in combinations(all_seg_keys, n_A_segs):
                shuffled_B_segs = tuple(s for s in all_seg_keys if s not in shuffled_A_segs)

                if (shuffled_A_segs, shuffled_B_segs) in seen or (shuffled_B_segs, shuffled_A_segs) in seen:
                    continue

                seen.add((shuffled_A_segs, shuffled_B_segs))
                seen.add((shuffled_B_segs, shuffled_A_segs))

                yield shuffled_A_segs, shuffled_B_segs

    background_scores = dict((g,
                              dict((k, {})
                                   for k in [OVERALL_SCORE] + states))
                             for g in groups)

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

                                          dict((groups[g_idx],
                                                dict((dataset,
                                                      all_segs[dataset]) for dataset in group_datasets))
                                               for g_idx, group_datasets in enumerate(shuffled_group_segmentations)),

                                          BIN_SIZE,
                                          states,
                                          chrom_lengths,

                                          max_min_window,
                                          to_smooth,

                                          use_mean_distance_matrix,
                                          compute_differential_from_closest_group
                                          )
                                         for combo_no, shuffled_group_segmentations in
                                         enumerate(shuffled_samples(n_perm)))):

        for g in groups:
            for score_type in combo_background_scores[g]:
                for s in combo_background_scores[g][score_type]:
                    if s not in background_scores[g][score_type]:
                        background_scores[g][score_type][s] = 0
                    background_scores[g][score_type][s] += combo_background_scores[g][score_type][s]
        i += 1
        echo('Main thread update:', i)
        gc.collect()

    if n_threads > 1:
        pool.close()

    for g in groups:
        for score_type in background_scores[g]:
            total_regions = sum(background_scores[g][score_type].itervalues())
            total_so_far = 0

            for s in sorted(background_scores[g][score_type], reverse=True):
                total_so_far += background_scores[g][score_type][s]
                background_scores[g][score_type][s] = total_so_far / float(total_regions)

    clear_metric_cache()

    return dict((g,
                 dict((score_type,
                       sorted(background_scores[g][score_type].items()))
                      for score_type in background_scores))
                for g in groups)


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


def process_shuffled_combo(combo_no,

                           shuffled_group_segmentations,

                           BIN_SIZE,
                           states,
                           chrom_lengths,
                           max_min_window,
                           to_smooth,

                           use_mean_distance_matrix,
                           compute_differential_from_closest_group):

    echo('Combo no:', combo_no)
    shuffled_metric = learn_metrics_for_groups(shuffled_group_segmentations, states, BIN_SIZE, use_mean_distance_matrix)

    # iterate over the bins in both segmentation groups

    groups = sorted(shuffled_group_segmentations)

    background_scores = dict((g, dict((k, {}) for k in [OVERALL_SCORE] + states)) for g in groups)

    random_chunks = generate_random_chunks(None,
                                           None,
                                           BIN_SIZE,
                                           N_RANDOM_CHUNKS=10,
                                           _chrom_lengths=chrom_lengths)

    for chrom in sorted(random_chunks):
        # echo(combo_no, chrom)
        chrom_segs = dict((g, dict((d,
                                    decompress_chrom(shuffled_group_segmentations[g][d][chrom])
                                    ))) for g in groups)

        for chunk_start, chunk_end in random_chunks[chrom]:
            # echo(combo_no, chunk_start, chunk_end)

            chunk_scores = get_overall_and_per_state_diff_score_C(chrom,
                                                                  dict((g,
                                                                        slice_segmentations(chrom_segs[g],
                                                                                            chunk_start,
                                                                                            chunk_end))
                                                                       for g in chrom_segs),
                                                                  shuffled_metric,
                                                                  states,
                                                                  None,
                                                                  True,
                                                                  max_min_window,
                                                                  background_chunk=True,
                                                                  use_posteriors=False,
                                                                  compute_differential_from_closest_group=compute_differential_from_closest_group)

            if to_smooth:
                smooth_dict(chunk_scores) #chunk_scores = dict((score_type, smooth(chunk_scores[score_type])) for score_type in chunk_scores)

            for g in background_scores:
                for score_type in background_scores[g]:
                    for s in chunk_scores[g][score_type]:

                        abs_score = abs(s)

                        if abs_score not in background_scores[g][score_type]:
                            background_scores[g][score_type][abs_score] = 0

                        background_scores[g][score_type][abs_score] += 1
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
                              posteriors_dir=None):

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
                                                                use_posteriors=(posteriors_dir is not None))

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
                             group_segmentations,
                             metric,
                             states,
                             use_mean_distance_matrix):

    MAX_N_SHUFFLED_SAMPLES = 100
    groups = sorted(group_segmentations)

    echo('Learning significance threshold by permuting segmentations at FDR:', fdr_threshold)

    background_scores = compute_background_scores_by_shuffling_segmentations(group_segmentations,

                                                                             states,

                                                                             metric,

                                                                             n_perm=MAX_N_SHUFFLED_SAMPLES,

                                                                             to_smooth=args.smooth,

                                                                             max_min_window=args.max_min_window,
                                                                             n_threads=args.n_threads,

                                                                             use_mean_distance_matrix=use_mean_distance_matrix,
                                                                             compute_differential_from_closest_group=args.compute_differential_from_closest_group)

    foreground_scores = compute_foreground_scores(group_segmentations,
                                                  states,
                                                  metric,

                                                  to_smooth=args.smooth,
                                                  max_min_window=args.max_min_window,
                                                  compute_differential_from_closest_group=args.compute_differential_from_closest_group
                                                  )

    background_model = dict((g, dict((k, {0.0: 1}) for k in [OVERALL_SCORE] + states)) for g in groups)
    background_threshold = dict((g, dict((k, None) for k in [OVERALL_SCORE] + states)) for g in groups)

    # determine FDR threshold
    for g in groups:
        for score_type in sorted(background_model[g]):

            all_scores = foreground_scores[g][score_type]

            best_score_cutoff = 0
            best_fdr = 1

            score_idx = 0
            bgr_idx = 0

            while score_idx < len(all_scores):

                # find how many background scores are greater than the current real score
                while bgr_idx < len(background_scores[g][score_type]) and background_scores[g][score_type][bgr_idx][0] < all_scores[score_idx][0]:
                    bgr_idx += 1

                if bgr_idx == len(background_scores[g][score_type]):
                    false_positives = 1
                else:
                    false_positives = background_scores[g][score_type][bgr_idx][1]

                all_positives = all_scores[score_idx][1]

                current_fdr = false_positives / all_positives

                if current_fdr < best_fdr:
                    best_fdr = current_fdr
                    best_score_cutoff = all_scores[score_idx][0]
                else:
                    current_fdr = best_fdr

                background_model[g][score_type][all_scores[score_idx][0]] = current_fdr

                if current_fdr <= fdr_threshold and background_threshold[g][score_type] is None:
                    background_threshold[g][score_type] = all_scores[score_idx][0]

                score_idx += 1

            echo(g, score_type,
                 'Score cutoff at FDR', fdr_threshold, ':', background_threshold[g][score_type],
                 'bgr model size:', len(background_model[g][score_type]))

            if background_threshold[g][score_type] is None:
                echo(g, 'No significant changes were found at FDR of', fdr_threshold, '. Best FDR=', best_fdr, 'cutoff=', best_score_cutoff)
                # background_threshold[score_type] = all_scores[-1][0] + 1

                # background_threshold[score_type] = best_score_cutoff

    return (dict((g,
                 dict((score_type,
                       sorted(background_model[g][score_type].iteritems()))
                      for score_type in background_model[g]))
                for g in groups),
            background_threshold)

SIGNAL = 'signal'
STATE_TRANSITIONS = 'state_transitions'


def decompress_chrom(chrom_data):
    dec_chrom = zlib.decompress(chrom_data).split('\n')
    return map(lambda s: parse(s, [int, int, str]), dec_chrom)


def process_chromosome(out_prefix,

                       chrom,
                       chrom_length,
                       group_segmentations,
                       metric,

                       states,
                       consistent_state_probs,
                       compute_per_state_scores,
                       max_min_window,
                       BIN_SIZE=None,
                       background_chunk=False,
                       use_posteriors=False,
                       report_per_group_state_probabilities=False,
                       constitutive_only=False,
                       compute_differential_from_closest_group=False):

    echo('Processing:', chrom, 'length:', chrom_length, 'bins')

    # chrom, signal, closest_state = get_overall_and_per_state_diff_score(chrom,
    chrom, signal, closest_state = get_overall_and_per_state_diff_score_C(chrom,
                                                                          chrom_length,

                                                                        group_segmentations,
                                                                        metric,

                                                                        states,
                                                                        consistent_state_probs,
                                                                        compute_per_state_scores,
                                                                        max_min_window,
                                                                        background_chunk,
                                                                        use_posteriors,
                                                                        report_per_group_state_probabilities,
                                                                        constitutive_only,
                                                                        compute_differential_from_closest_group)

    echo('Writing down pickle files:', chrom)
    gc.collect()

    pickle_fnames = {}
    st_fnames = {}

    for group in signal:

        if args.smooth:
            echo('Smoothing:', group, chrom)
            smooth_dict(signal[group])

        st_fname = out_prefix + '.' + group + '.' + chrom + '.' + STATE_TRANSITIONS + '.pickle'
        st_fnames[group] = st_fname

        with open(st_fname, 'w') as out_f:
            pickle.dump(closest_state[group], out_f, pickle.HIGHEST_PROTOCOL)

        pickle_fnames[group] = {}
        for score_type in signal[group]:
            pickle_fname = out_prefix + '.' + group + '.' + chrom + '.' + score_type + '.pickle'
            pickle_fnames[group][score_type] = pickle_fname

            with open(pickle_fname, 'w') as out_f:
                pickle.dump(signal[group][score_type], out_f, pickle.HIGHEST_PROTOCOL)

    gc.collect()
    echo('Done:', chrom)

    return pickle_fnames, st_fnames, chrom


def merge_output_files(out_prefix,
                       group,
                       score_type,
                       score_fnames,
                       st_fnames,
                       BIN_SIZE,
                       score_background_model,
                       threshold,
                       skip_bed_files=False):

    echo('Merging:', group, score_type)
    wig_fname = out_prefix + '.' + group + '.' + re.sub(r'\W+', '_', score_type) + '.wig.gz'
    bed_fname = out_prefix + '.' + group + '.' + re.sub(r'\W+', '_', score_type) + '.summary.bed.gz'

    for fname in [wig_fname, bed_fname]:
        if os.path.isfile(fname):
            os.unlink(fname)

    echo('Reading pickle files:', group, score_type)

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
    sorted_scores_indexes = dict((c, 0) for c in signal)
    sorted_scores = {}
    for chrom in signal.keys():
        sorted_scores[chrom] = sorted([(score, bin_no, state_transition)
                                            for bin_no, (score, state_transition) in enumerate(izip(signal[chrom],
                                                                                                    state_transitions[chrom]))],
                                      key=lambda (score, bin_no, state_transition): (-score, bin_no))

        del signal[chrom]
        del state_transitions[chrom]
    gc.collect()
    # sorted_scores = sorted([(score, chrom, bin_no, state_transition)
    #                                 for chrom in signal
    #                                     for bin_no, (score, state_transition) in enumerate(izip(signal[chrom],
    #                                                                                             state_transitions[chrom]))],
    #                        key=lambda (score, chrom, bin_no, state_transition): (-score, chrom, bin_no))

    echo('Writing sorted summary files for', score_type, ':', bed_fname)

    score_cache = {}
    def get_score(chrom, sorted_scores, sorted_scores_indexes):
        idx = sorted_scores_indexes[chrom]

        if idx == len(sorted_scores[chrom]):
            return None

        sorted_scores_indexes[chrom] += 1
        return sorted_scores[chrom][idx], chrom

    chromosomes = sorted(sorted_scores)

    c_scores = dict((chrom, get_score(chrom, sorted_scores, sorted_scores_indexes)) for chrom in chromosomes)

    with open_file(bed_fname, 'w') as out_f:
        while len(c_scores) > 0:
            (score, bin_no, closest_state), chrom = max(c_scores.itervalues())

        # for score, chrom, bin_no, (closest_state_A, closest_state_B, a_states, b_states) in sorted_scores:
            out_f.write('%s\t%d\t%d\t%lf\t%s\t%lf\t%s\n' %
                        (chrom,
                         bin_no * BIN_SIZE,
                         (bin_no + 1) * BIN_SIZE,
                         score,

                         closest_state,

                         get_FDR(abs(score), score_background_model, score_cache),
                         '1' if threshold is not None and abs(score) >= threshold else '0'))

            new_score = get_score(chrom, sorted_scores, sorted_scores_indexes)

            if new_score is None:
                del c_scores[chrom]
            else:
                c_scores[chrom] = new_score

    gc.collect()

    echo('Done:', score_type)
    return True


def _compute_average_metric(metric, states):
    average_metric = dict((s1, dict((s2, 0) for s2 in states)) for s1 in states)

    total_groups = len(metric)

    for g in metric:
        group_average_metric = dict((s1, dict((s2, 0) for s2 in states)) for s1 in states)
        m = metric[g]
        for d in m:
            for s1 in m[d]:
                for s2 in m[d][s1]:
                    group_average_metric[s1][s2] += m[d][s1][s2]

        for s1 in group_average_metric:
            for s2 in group_average_metric[s1]:
                average_metric[s1][s2] += group_average_metric[s1][s2] / len(metric[g])

    average_metric = dict((s1, dict((s2, average_metric[s1][s2] / total_groups) for s2 in states)) for s1 in states)

    return dict((g, dict((d, average_metric) for d in metric[g])) for g in metric)


def average_distribution(distributions, states):
    avg_dist = dict((s, 0) for s in states)
    n = 0
    for d in distributions:
        n += 1
        for s in d:
            avg_dist[s] += d[s]

    for s in avg_dist:
        avg_dist[s] /= n

    return avg_dist


def _print_consistent_state_transition_scores(states, metric_A, metric_B):
    to_print = '\nRaw' + '\n'
    to_print += '\t'.join(['State'] + states) + '\n'
    for s1 in states:
        a_probs = get_prob_vector([s1 for _ in metric_A], metric_A.keys(), metric_A, states)
        to_print += '\t'.join([s1] + [str(symmetric_KL_divergence(
                                            a_probs,
                                            average_distribution(
                                                (get_prob_vector([s2 for _ in metric_B[g]],
                                                                 metric_B[g].keys(),
                                                                 metric_B[g],
                                                                 states) for g in metric_B),
                                                states)))
                                for s2 in states]) + '\n'

    echo(to_print)


def _generate_consistent_state_probabilities(states, metric):

    consistent_state_probs = {}

    for g in metric:
        consistent_state_probs[g] = {}
        for s in states:
            consistent_state_probs[g][s] = get_prob_vector([s for _ in metric[g]],
                                                                       metric[g].keys(),
                                                                       metric[g],
                                                                       states)

    return consistent_state_probs


def decompress_group(group_segs):
    dec = {}
    for d in group_segs:
        dec[d] = {}
        for chrom in group_segs[d]:
            dec_chrom = zlib.decompress(group_segs[d][chrom]).split('\n')
            dec[d][chrom] = map(lambda s: parse(s, [int, int, str]), dec_chrom)

    return dec


def learn_metrics_for_groups(group_segmentations, states, BIN_SIZE, use_mean_distance_matrix):
    metric = {}

    for g in group_segmentations:
        if len(group_segmentations[g]) > 1:
            echo('Learning metrics for:', g)
            metric[g] = learn_metric_from_all_replicates(decompress_group(group_segmentations[g]),
                                                         states,
                                                         BIN_SIZE,
                                                         None)

    average_metric = _compute_average_metric(metric, states)

    if use_mean_distance_matrix:
        metric = average_metric

    else:

        single_avg_metric = average_metric.values()[0].values()[0]

        for g in group_segmentations:
            if len(group_segmentations[g]) == 1:
                echo('Using average metric for:', g)
                fname = group_segmentations[g].keys()[0]
                metric[g] = {fname: single_avg_metric}

    return metric

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-g', nargs='+', help='files with paths to segmentation files for all groups')

    parser.add_argument('--bin-size', type=int, help='Bin size, default: %(default)s', default=200)
    parser.add_argument('-f', '--fdr-threshold', type=float, help='FDR threshold: %(default)s', default=0.05)

    parser.add_argument('-o', '--output', help='output prefix')
    # parser.add_argument('--background-scores', help='bacground_scores.pickle')

    # parser.add_argument('--chrom-hmm-binarized', help='path to ChromHMM binarized files to compute background scores')
    # parser.add_argument('--chrom-hmm-model-path', help='path to the ChromHMM model to use for shuffled marks')

    parser.add_argument('--posteriors-dir', dest='posteriors_dir', help='posteriors directory from ChromHMM (scores will be computed by'
                                                                        'taking into account the posteriors')

    parser.add_argument('-s', '--smooth', action='store_true', help='smooth signal')
    parser.add_argument('-p', '--n-threads', type=int, help='number of threads to use', default=1)

    parser.add_argument('--chroms', nargs='+', help='use only these chromosomes', default=None)

    parser.add_argument('--max-min-window', type=int, help='pick the maximum distance between the two most similar bins within this window', default=0)

    parser.add_argument('--metrics', dest='metric_fname', help='Pickle file with stored metrics', default=None)

    # parser.add_argument('--per-state-scores', action='store_true', help='compute per state scores', default=True)
    parser.add_argument('--constitutive-only', action='store_true', help='constitutive scores only', default=False)
    parser.add_argument('--report-per-group-state-probabilities', action='store_true', help='output the probabilities for each state in each group per state scores', default=False)
    parser.add_argument('--skip-bed-files', action='store_true', help='output only wiggle files', default=False)

    parser.add_argument('--use-mean-distance-matrix', action='store_true',
                        help='Use one distance matrix for all samples computed as the average distance between each '
                             'pair of states',
                        default=False)

    parser.add_argument('--compute-differential-from-closest-group',
                        dest='compute_differential_from_closest_group',
                        action='store_true',
                        help='Use the closest group to compute differential scores [%(default)s]',
                        default=False)

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)

    out_fname = args.output

    open_log(out_fname + '.log')
    BIN_SIZE = args.bin_size
    echo('cmd:', ' '.join(sys.argv))
    for arg in sorted(vars(args)):
        echo(arg, '=', getattr(args, arg))

    posteriors_dir = args.posteriors_dir
    use_mean_distance_matrix = args.use_mean_distance_matrix

    max_min_window = args.max_min_window

    # compute_per_state_scores = args.per_state_scores

    get_seg_fnames = lambda fname: [os.path.join(os.path.split(fname)[0], l.strip()) for l in open_file(fname).readlines()]

    group_segmentations = {}

    states = None

    chrom_lengths = {}
    for g in args.g:
        echo('Group:', g)
        # segs, states = read_segmentations(get_seg_fnames(g))
        segs, states, cur_chrom_lengths = read_segmentations_compressed(get_seg_fnames(g))

        group_segmentations[os.path.split(g)[1].replace('.seg_files', '').replace('.segmentations', '')] = segs

        for chrom in cur_chrom_lengths:
            if chrom not in chrom_lengths:
                chrom_lengths[chrom] = cur_chrom_lengths[chrom] / BIN_SIZE
            else:
                chrom_lengths[chrom] = min(chrom_lengths[chrom], cur_chrom_lengths[chrom] / BIN_SIZE)

    echo('Chromosome Lengths:', chrom_lengths)
    if args.chroms is not None:
        chromosomes = args.chroms
    else:
        chromosomes = sorted(reduce(operator.and_, [set(s.keys())
                                             for segs in group_segmentations.values()
                                             for s in segs.values()]),
                         key=lambda c: int(c.replace('chr','')) if re.search(r'^chr(\d+)$', c) else 100)

    echo('Chromosomes:', chromosomes)
    # chromosomes = ['chr22']
    # print chromosomes
    group_segmentations = dict((g, filter_chroms(group_segmentations[g], chromosomes)) for g in group_segmentations)

    # print 'A'
    # metric_A = learn_metric_from_all_replicates_test(group_A_segmentations, states, BIN_SIZE, posteriors_dir)
    # print 'B'
    # metric_B = learn_metric_from_all_replicates_test(group_B_segmentations, states, BIN_SIZE, posteriors_dir)
    # exit(1)


    if args.metric_fname is not None:
        echo('Loading metrics from:', args.metric_fname)
        with open(args.metric_fname) as in_f:
            metric = pickle.load(in_f)

    else:
        metric = learn_metrics_for_groups(group_segmentations, states, BIN_SIZE, use_mean_distance_matrix)

        with open(out_fname + '.metrics.pickle', 'w') as m_out:
            echo('Storing metrics in:', out_fname + '.metrics.pickle')
            pickle.dump(metric, m_out, protocol=pickle.HIGHEST_PROTOCOL)

    echo('Metrics')
    for g in metric:
        echo('Group:', g)
        print_metric(metric[g])

    score_types = [OVERALL_SCORE] + states # if compute_per_state_scores else [])

    # _print_consistent_state_transition_scores(states, metric_A, metric_B)

    fdr_threshold = args.fdr_threshold

    if 0 < fdr_threshold < 1:
        background_model, background_threshold = compute_background_model(args,
                                                                          fdr_threshold,
                                                                          group_segmentations,
                                                                          metric,
                                                                          states,
                                                                          use_mean_distance_matrix)
        gc.collect()
    else:
        echo('Skipping null model computations')
        background_model = dict((g, dict((s, [(0, 0)]) for s in score_types)) for g in list(group_segmentations) + [CONSTITUTIVE])
        background_threshold = dict((g, dict((s, 0) for s in score_types)) for g in list(group_segmentations) + [CONSTITUTIVE])

    real_scores = {}
    real_state_transitions = {}
    echo('Computing real scores')

    consistent_state_probs = _generate_consistent_state_probabilities(states, metric)

    if args.n_threads > 1:
        pool = Pool(args.n_threads)
        _map = pool.imap_unordered

    score_fnames = {}
    st_fnames = {}
    for pickle_fnames, st_fname, chrom in _map(worker,
                                                 ((process_chromosome,

                                                   out_fname,
                                                   _chrom,
                                                   chrom_lengths[_chrom],

                                                   dict((g,
                                                         dict((d, group_segmentations[g][d][_chrom])
                                                        for d in group_segmentations[g])) for g in group_segmentations),
                                                   # dict((g,
                                                   #       dict((d, chrom_segmentation_to_list(group_segmentations[g][d][_chrom], BIN_SIZE)
                                                   #               if posteriors_dir is None else
                                                   #             read_posteriors_as_list(get_chrom_posteriors_fname(posteriors_dir, d, _chrom)))
                                                   #      for d in group_segmentations[g])) for g in group_segmentations),

                                                   metric,

                                                   states,
                                                   consistent_state_probs,
                                                   True,
                                                   max_min_window,
                                                   BIN_SIZE,
                                                   False,
                                                   (posteriors_dir is not None),
                                                   args.report_per_group_state_probabilities,
                                                   args.constitutive_only,
                                                   args.compute_differential_from_closest_group) for _chrom in chromosomes)):

        for group in pickle_fnames:

            if group not in st_fnames:
                st_fnames[group] = {}

            st_fnames[group][chrom] = st_fname[group]

            if group not in score_fnames:
                score_fnames[group] = {}

            for score_type in pickle_fnames[group]:
                if score_type not in score_fnames[group]:
                    score_fnames[group][score_type] = {}

                score_fnames[group][score_type][chrom] = pickle_fnames[group][score_type]

    # print score_fnames
    # print st_fnames

    all(_map(worker, ((merge_output_files,
                       out_fname,
                       group,
                       score_type,
                       score_fnames[group][score_type],
                       st_fnames[group],
                       BIN_SIZE,
                       background_model[group].get(score_type, [(0, 0)]),
                       background_threshold[group].get(score_type, None),
                       args.skip_bed_files
                       ) for group in sorted(score_fnames)
                            for score_type in sorted(score_fnames[group]))))

    # delete pickled state transition files
    for group in st_fnames:
        for fname in st_fnames[group].values():
            os.unlink(fname)

    if args.n_threads > 1:
        pool.close()

    echo('Done')
    close_log()
