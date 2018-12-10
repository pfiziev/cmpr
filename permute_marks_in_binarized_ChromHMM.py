import argparse
import gzip
from itertools import izip
import os
import shutil
import sys, random
import tempfile
import operator
import re
import gc
from ChromDelta import get_overall_and_per_state_diff_score, worker
from metric import *
from utils import echo, read_segmentations, filter_chroms, smooth
from subprocess import call
import zlib
import itertools
from multiprocessing import Pool

CHROM_HMM_EXEC = '/Users/pfiziev/software/ChromHMM/ChromHMM.jar'
# CHROM_HMM_EXEC = '/u/home/p/pfiziev/software/ChromHMM/ChromHMMp.jar'


def read_binarized(fname, random_chunks):
    echo('Reading:', fname)
    with open_file(fname) as in_f:
        ct, chrom = in_f.readline().strip().split()

        if chrom not in random_chunks:
            return None, None, None, None

        marks = in_f.readline().strip().split()
        c_data = {}

        chunks = random_chunks[chrom]
        chunk_no = 0
        chunk_start, chunk_end = chunks[chunk_no]
        chunk_id = chrom + '_' + str(chunk_no)
        c_data[chunk_id] = dict((m, []) for m in marks)

        for bin_no, l in enumerate(in_f):

            if bin_no == chunk_end:
                chunk_no += 1
                if chunk_no == len(chunks):
                    break

                chunk_id = chrom + '_' + str(chunk_no)
                c_data[chunk_id] = dict((m, []) for m in marks)
                chunk_start, chunk_end = chunks[chunk_no]

            if bin_no < chunk_start:
                continue

            binarized = l.strip().split()
            for b, m in zip(binarized, marks):
                c_data[chunk_id][m].append(b)

        # convert the binary data array to a string to save some memory
        for chunk_id in c_data:
            for m in c_data[chunk_id]:
                c_data[chunk_id][m] = ''.join(c_data[chunk_id][m])

    return marks, ct, chrom, c_data


def str_to_bits(string_array):
    bits = [0] * (len(string_array) / 8 + (1 if len(string_array) % 8 != 0 else 0))
    c_bit_array_idx = -1
    for i in xrange(len(string_array)):
        bit_pos = i % 8
        if bit_pos == 0:
            c_bit_array_idx += 1
        if string_array[i] not in ['0', '1']:
            echo('ERROR: ChromHMM binarized data should consist of only 0 and 1s:', i, ', char=', string_array[i])
            exit(1)
        bits[c_bit_array_idx] |= int(string_array[i]) << bit_pos
    return ''.join(map(chr, bits))


def bits_to_string(bit_string, original_length):
    string_array = [''] * original_length
    string_idx = 0
    for bit in bit_string:
        for i in xrange(8):
            string_array[string_idx] = str(int(bool(ord(bit) & (1 << i))))
            string_idx += 1
            if string_idx == original_length:
                break
    return ''.join(string_array)


def compress(array):
    return zlib.compress(''.join(array))
    # return zlib.compress(str_to_bits(array))


def decompress(compressed, original_length):
    return zlib.decompress(compressed)
    # return bits_to_string(zlib.decompress(compressed), original_length)


def read_binarized_in_full(fname, chrom_lengths):
    echo('Reading:', fname)
    with open_file(fname) as in_f:
        ct, chrom = in_f.readline().strip().split()

        # if chrom not in chrom_lengths:
        #     echo('Skipping:', chrom)
        #     return None, None, None, None
        #
        marks = in_f.readline().strip().split()
        data = dict((m, [''] * chrom_lengths[chrom]) for m in marks)

        for bin_no, l in enumerate(in_f):
            binarized = l.strip().split()
            for b, m in zip(binarized, marks):
                data[m][bin_no] = b

        for m in data:
            data[m] = compress(data[m])

    return marks, ct, chrom, data


def clean_up_directory(dirname):
    for f in os.listdir(dirname):
        os.remove(os.path.join(dirname, f))
#
# def _compute_background_scores_by_shuffling_marks(bin_directory,
#                                                  seg_fnames,
#                                                  ChromHMM_model_path,
#                                                  real_chrom_lengths,
#                                                  n_group_A,
#                                                  BIN_SIZE,
#                                                  states,
#                                                  n_perms=100,
#                                                  compute_per_state_scores=False,
#                                                  max_min_window=0,
#                                                  to_smooth=False,
#                                                  use_posteriors=False):
#
#     ct_chrom_fnames = {}
#     celltypes = set(get_celltype(sf) for sf in seg_fnames)
#
#     out_dir = tempfile.mkdtemp()
#     echo('temp dir:', out_dir)
#
#     data = {}
#     marks = None
#     background_model = new_score_dict(states, compute_per_state_scores)
#
#
#     # chunk_chrom_lengths = {}
#     for fname in sorted(os.listdir(bin_directory)):
#         if re.search(r'_binary\.txt(\.gz)?', fname) and any(ct in fname for ct in celltypes):
#             fname = os.path.join(bin_directory, fname)
#             c_marks, ct, chrom, c_data = read_binarized_in_full(fname, real_chrom_lengths)
#
#             if marks is None:
#                 marks = c_marks
#
#             if ct not in data:
#                 data[ct] = {}
#                 ct_chrom_fnames[ct] = {}
#
#             data[ct][chrom] = c_data
#             ct_chrom_fnames[ct][chrom] = fname
#
#     mark_celltype_order = dict((m, list(sorted(data.keys()))) for m in marks)
#
#     seen = set()
#     # seen.add(tuple([d for m in mark_celltype_order for d in mark_celltype_order[m]]))
#     for perm_idx in xrange(n_perms):
#         if perm_idx >= 0:
#             key = None
#             while key is None or key in seen:
#                 for m in mark_celltype_order:
#                     random.shuffle(mark_celltype_order[m])
#                 key = tuple([d for m in sorted(mark_celltype_order) for d in mark_celltype_order[m]])
#             seen.add(key)
#         echo('Combo:', perm_idx)
#         echo('Writing shuffled binarized files')
#
#         random_chunks = generate_random_chunks(None, None, BIN_SIZE, N_RANDOM_CHUNKS=100,
#                                                _chrom_lengths=real_chrom_lengths)
#
#         for chrom in sorted(random_chunks):
#             decompressed_data = dict((d,
#                                       dict((m,
#                                             zlib.decompress(data[d][chrom][m]))
#                                            for m in marks))
#                                      for d in data)
#
#             for chunk_no, (chunk_start, chunk_end) in enumerate(random_chunks[chrom]):
#                 chunk_length = chunk_end - chunk_start
#                 chunk_chrom = chrom + '_' + str(chunk_no)
#                 for ct_idx, ct in enumerate(sorted(celltypes)):
#                     _fname = ct_chrom_fnames[ct][chrom].replace('_' + chrom + '_',
#                                                                 '_' + chunk_chrom + '_')
#                     out_fname = re.sub(r'\.gz$', '', os.path.join(out_dir, os.path.split(_fname)[1]))
#
#                     # echo('marks:', marks)
#                     # echo('ct order:', [mark_celltype_order[m][ct_idx] for m in marks])
#
#                     # echo('Writing to:', out_fname)
#                     # chunk_data = [zlib.decompress(data[mark_celltype_order[m][ct_idx]][chrom][m])[chunk_start: chunk_end]
#                     #               for m in marks]
#                     chunk_data = [decompressed_data[mark_celltype_order[m][ct_idx]][m][chunk_start: chunk_end]
#                                   for m in marks]
#
#                     with open(out_fname, 'w') as out_f:
#
#                         out_f.write('\t'.join([ct, chunk_chrom]) + '\n')
#                         out_f.write('\t'.join(marks) + '\n')
#
#                         for bin_idx in xrange(chunk_length):
#                             out_f.write(
#                                 # '\t'.join([data[mark_celltype_order[m][ct_idx]][chrom][m][bin_idx] for m in marks]) + '\n')
#                                 '\t'.join([chunk_data[m_idx][bin_idx] for m_idx in xrange(len(marks))]) + '\n')
#
#         echo('Running ChromHMM')
#         if use_posteriors:
#             call(['java', '-mx2000M', '-jar', CHROM_HMM_EXEC, 'MakeSegmentation', '-printposterior', ChromHMM_model_path, out_dir, out_dir])
#         else:
#             call(['java', '-mx2000M', '-jar', CHROM_HMM_EXEC, 'MakeSegmentation', ChromHMM_model_path, out_dir, out_dir])
#
#         process_shuffled_segmentations(n_group_A,
#                                        out_dir,
#                                        background_model,
#                                        BIN_SIZE,
#                                        compute_per_state_scores,
#                                        max_min_window,
#                                        to_smooth,
#                                        use_posteriors)
#
#         clean_up_directory(out_dir)
#
#     echo('Purging temp directory:', out_dir)
#     shutil.rmtree(out_dir)
#
#     # convert the score counts to a CDF
#     for score_type in background_model:
#         total_regions = sum(background_model[score_type].itervalues())
#         total_so_far = 0
#
#         for s in sorted(background_model[score_type], reverse=True):
#             total_so_far += background_model[score_type][s]
#             background_model[score_type][s] = total_so_far / float(total_regions)
#
#     return dict((score_type, sorted(background_model[score_type].items())) for score_type in background_model)


def compute_background_scores_by_shuffling_marks( bin_directory,
                                                  seg_fnames,
                                                  ChromHMM_model_path,
                                                  real_chrom_lengths,

                                                  n_group_A,
                                                  BIN_SIZE,
                                                  states,
                                                  n_perms=100,
                                                  compute_per_state_scores=False,
                                                  max_min_window=0,
                                                  to_smooth=False,
                                                  use_posteriors=False,
                                                  n_threads=1,
                                                  use_mean_distance_matrix=False,
                                                  keep_max_scores_per_bin=False):

    ct_chrom_fnames = {}
    celltypes = set(get_celltype(sf) for sf in seg_fnames)

    # out_dir = tempfile.mkdtemp()
    # echo('temp dir:', out_dir)

    data = {}
    marks = None
    background_model = new_score_dict(states, compute_per_state_scores)


    # chunk_chrom_lengths = {}
    for fname in sorted(os.listdir(bin_directory)):
        if re.search(r'_binary\.txt(\.gz)?', fname) and \
                any((ct + '_') in fname for ct in celltypes) and \
                any(('_' + chrom + '_') in fname for chrom in real_chrom_lengths):
            fname = os.path.join(bin_directory, fname)
            c_marks, ct, chrom, c_data = read_binarized_in_full(fname, real_chrom_lengths)

            if marks is None:
                marks = c_marks

            if ct not in data:
                data[ct] = {}
                ct_chrom_fnames[ct] = {}

            data[ct][chrom] = c_data
            ct_chrom_fnames[ct][chrom] = fname

    gc.collect()

    mark_celltype_order = dict((m, list(sorted(data.keys()))) for m in marks)

    permutations = set()
    # seen.add(tuple([d for m in mark_celltype_order for d in mark_celltype_order[m]]))

    for perm_idx in xrange(n_perms):
        key = None
        while key is None or key in permutations:
            for m in mark_celltype_order:
                random.shuffle(mark_celltype_order[m])

            # key = tuple([d for m in sorted(mark_celltype_order) for d in mark_celltype_order[m]])

            key = tuple(tuple(mark_celltype_order[m]) for m in marks)

        permutations.add(key)

    if n_threads == 1:
        _map = itertools.imap
    else:
        pool = Pool(processes=n_threads)
        _map = pool.imap_unordered

    for i, combo_background_scores in enumerate(_map(worker,
                                                    ((process_a_set_of_marks_permutation,
                                                      ps_idx,
                                                      permutations_subset,
                                                      data,

                                                      n_group_A,

                                                      celltypes,
                                                      marks,
                                                      ct_chrom_fnames,

                                                      real_chrom_lengths,
                                                      ChromHMM_model_path,

                                                      BIN_SIZE,
                                                      states,
                                                      compute_per_state_scores,
                                                      max_min_window,
                                                      to_smooth,
                                                      use_posteriors,
                                                      use_mean_distance_matrix,
                                                      keep_max_scores_per_bin)

                                                    for ps_idx, permutations_subset
                                                            in enumerate(chunks(permutations,
                                                                                max(1, n_perms / n_threads)))))):

        for score_type in combo_background_scores:
            for s in combo_background_scores[score_type]:
                if s not in background_model[score_type]:
                    background_model[score_type][s] = 0

                background_model[score_type][s] += combo_background_scores[score_type][s]

        echo('Main thread combo update:', i)
        #gc.collect()

    if n_threads > 1:
        pool.close()

    # convert the score counts to a CDF
    for score_type in background_model:
        total_regions = sum(background_model[score_type].itervalues())
        total_so_far = 0

        for s in sorted(background_model[score_type], reverse=True):
            total_so_far += background_model[score_type][s]
            background_model[score_type][s] = total_so_far / float(total_regions)

    return dict((score_type, sorted(background_model[score_type].items())) for score_type in background_model)


def process_a_set_of_marks_permutation(ps_idx,
                                       permutations,
                                       data,

                                       n_group_A,
                                       celltypes,
                                       marks,

                                       ct_chrom_fnames,

                                       real_chrom_lengths,

                                       ChromHMM_model_path,
                                       BIN_SIZE,

                                       states,

                                       compute_per_state_scores,
                                       max_min_window,
                                       to_smooth,
                                       use_posteriors,
                                       use_mean_distance_matrix,
                                       keep_max_scores_per_bin):

    background_model = new_score_dict(states, compute_per_state_scores)

    for perm_no, perm in enumerate(permutations):
        echo('Subset:', ps_idx, 'Combo:', perm_no)
        combo_background_scores = process_marks_permutation(perm,
                                                            data,

                                                            n_group_A,

                                                            celltypes,
                                                            marks,
                                                            ct_chrom_fnames,

                                                            real_chrom_lengths,
                                                            ChromHMM_model_path,

                                                            BIN_SIZE,
                                                            states,
                                                            compute_per_state_scores,
                                                            max_min_window,
                                                            to_smooth,
                                                            use_posteriors,
                                                            use_mean_distance_matrix,
                                                            keep_max_scores_per_bin)
        for score_type in combo_background_scores:
            for s in combo_background_scores[score_type]:
                if s not in background_model[score_type]:
                    background_model[score_type][s] = 0

                background_model[score_type][s] += combo_background_scores[score_type][s]

        gc.collect()

    return background_model


def process_marks_permutation(perm,
                              data,

                              n_group_A,

                              celltypes,
                              marks,
                              ct_chrom_fnames,

                              real_chrom_lengths,
                              ChromHMM_model_path,

                              BIN_SIZE,
                              states,
                              compute_per_state_scores,
                              max_min_window,
                              to_smooth,
                              use_posteriors,
                              use_mean_distance_matrix,
                              keep_max_scores_per_bin):

    # echo('Combo:', perm_no)
    echo(perm)

    out_dir = tempfile.mkdtemp()
    echo('Writing shuffled binarized files to:', out_dir)

    background_model = new_score_dict(states, compute_per_state_scores)

    random_chunks = generate_random_chunks(None, None, BIN_SIZE, N_RANDOM_CHUNKS=100,
                                           _chrom_lengths=real_chrom_lengths)

    for chrom in sorted(random_chunks):
        decompressed_data = dict((d,
                                  dict((m,
                                        decompress(data[d][chrom][m], real_chrom_lengths[chrom]))
                                       for m in marks))
                                 for d in data)

        for chunk_no, (chunk_start, chunk_end) in enumerate(random_chunks[chrom]):
            chunk_length = chunk_end - chunk_start
            chunk_chrom = chrom + '_' + str(chunk_no)
            for ct_idx, ct in enumerate(sorted(celltypes)):
                _fname = ct_chrom_fnames[ct][chrom].replace('_' + chrom + '_',
                                                            '_' + chunk_chrom + '_')
                out_fname = re.sub(r'\.gz$', '', os.path.join(out_dir, os.path.split(_fname)[1]))

                chunk_data = [decompressed_data[perm[m_idx][ct_idx]][m][chunk_start: chunk_end]
                                    for m_idx, m in enumerate(marks)]

                with open(out_fname, 'w') as out_f:

                    out_f.write('\t'.join([ct, chunk_chrom]) + '\n')
                    out_f.write('\t'.join(marks) + '\n')

                    for bin_idx in xrange(chunk_length):
                        out_f.write('\t'.join([chunk_data[m_idx][bin_idx] for m_idx in xrange(len(marks))]) + '\n')

    echo('Running ChromHMM')
    if use_posteriors:
        call(['java', '-mx2000M', '-jar', CHROM_HMM_EXEC, 'MakeSegmentation', '-printposterior', ChromHMM_model_path, out_dir, out_dir])
    else:
        call(['java', '-mx2000M', '-jar', CHROM_HMM_EXEC, 'MakeSegmentation', ChromHMM_model_path, out_dir, out_dir])

    process_shuffled_segmentations(states,
                                   n_group_A,
                                   out_dir,
                                   background_model,
                                   BIN_SIZE,
                                   compute_per_state_scores,
                                   max_min_window,
                                   to_smooth,
                                   use_posteriors,
                                   use_mean_distance_matrix,
                                   keep_max_scores_per_bin)

    echo('Purging temp directory:', out_dir)
    shutil.rmtree(out_dir)

    return background_model

# echo('Purging temp directory:', out_dir)
#     shutil.rmtree(out_dir)
#
#     # convert the score counts to a CDF
#     for score_type in background_model:
#         total_regions = sum(background_model[score_type].itervalues())
#         total_so_far = 0
#
#         for s in sorted(background_model[score_type], reverse=True):
#             total_so_far += background_model[score_type][s]
#             background_model[score_type][s] = total_so_far / float(total_regions)
#
#     return dict((score_type, sorted(background_model[score_type].items())) for score_type in background_model)


def process_shuffled_segmentations(states,
                                   n_group_A,
                                   out_dir,
                                   background_model,
                                   BIN_SIZE,
                                   compute_per_state_scores,
                                   max_min_window,
                                   to_smooth=False,
                                   use_posteriors=False,
                                   use_mean_distance_matrix=False,
                                   keep_max_scores_per_bin=False):

    seg_fnames = sorted([os.path.join(out_dir, f) for f in os.listdir(out_dir) if f.endswith('_segments.bed')])

    group_A_seg_fnames = seg_fnames[:n_group_A]
    group_B_seg_fnames = seg_fnames[n_group_A:]

    group_A_segmentations, _ = read_segmentations(group_A_seg_fnames)
    group_B_segmentations, _ = read_segmentations(group_B_seg_fnames)

    chromosomes = sorted(reduce(operator.and_, [set(s.keys())
                                             for segs in [group_A_segmentations, group_B_segmentations]
                                             for s in segs.values()]),
                         key=lambda c: int(c.replace('chr','')) if re.search(r'^chr(\d+)$', c) else 100)
    # print chromosomes
    group_A_segmentations = filter_chroms(group_A_segmentations, chromosomes)
    group_B_segmentations = filter_chroms(group_B_segmentations, chromosomes)

    posteriors_dir = os.path.join(out_dir, 'POSTERIOR') if use_posteriors else None

    metric_A = learn_metric_from_all_replicates(group_A_segmentations, states, BIN_SIZE, posteriors_dir=posteriors_dir, clear_cache=True)
    # print_metric(metric_A)

    metric_B = learn_metric_from_all_replicates(group_B_segmentations, states, BIN_SIZE, posteriors_dir=posteriors_dir, clear_cache=True)
    # print_metric(metric_B)

    if use_mean_distance_matrix:
        metric_A, metric_B = compute_average_metric(metric_A, metric_B, states)

    # print_average_metric(states, metric_A, metric_B)
    N_RANDOM_CHUNKS = 10
    for chrom in random.sample(chromosomes, N_RANDOM_CHUNKS):
        if use_posteriors:
            chrom_segmentations_A = dict((d, read_posteriors_as_list(get_chrom_posteriors_fname(posteriors_dir, d, chrom))) for d in group_A_seg_fnames)
            chrom_segmentations_B = dict((d, read_posteriors_as_list(get_chrom_posteriors_fname(posteriors_dir, d, chrom))) for d in group_B_seg_fnames)
        else:
            chrom_segmentations_A = dict((d, chrom_segmentation_to_list(group_A_segmentations[d][chrom], BIN_SIZE)) for d in group_A_seg_fnames)
            chrom_segmentations_B = dict((d, chrom_segmentation_to_list(group_B_segmentations[d][chrom], BIN_SIZE)) for d in group_B_seg_fnames)


        chunk_scores = get_overall_and_per_state_diff_score(chrom,
                                                            chrom_segmentations_A,
                                                            metric_A,
                                                            chrom_segmentations_B,
                                                            metric_B,
                                                            states,
                                                            None,
                                                            compute_per_state_scores,
                                                            max_min_window,
                                                            background_chunk=True,
                                                            use_posteriors=use_posteriors,
                                                            keep_max_scores_per_bin=keep_max_scores_per_bin)
        echo('Keys:', sorted(chunk_scores.keys()))
        if to_smooth:
            smooth_dict(chunk_scores)   # chunk_scores = dict((score_type, smooth(chunk_scores[score_type])) for score_type in chunk_scores)

        for score_type in background_model:
            for s in chunk_scores[score_type]:
                abs_score = abs(s)

                if abs_score not in background_model[score_type]:
                    background_model[score_type][abs_score] = 0

                background_model[score_type][abs_score] += 1

#
# if __name__ == '__main__':
#
#     parser = argparse.ArgumentParser()
#
#     parser.add_argument('-i', nargs='+', help='directory with binarized chromosome files')
#     parser.add_argument('-s', nargs='+', help='original chromhmm segmentations')
#
#     parser.add_argument('-n', help='n permutations', type=int)
#     parser.add_argument('-a', help='number of samples in group A', type=int)
#     parser.add_argument('--bin-size', default=200, help='bin size', type=int)
#     parser.add_argument('-m', help='ChromHMM model')
#
#     args = parser.parse_args()
#
#     if len(sys.argv) == 1:
#         parser.print_help()
#         exit(1)
#
#     n_perms = args.n
#     bin_directory = args.i
#     seg_fnames = args.s
#     ChromHMM_model_path = args.m
#     n_group_A = args.a
#     BIN_SIZE = args.bin_size
#
#     # compute_background_scores_by_shuffling_marks(bin_directory, seg_fnames, ChromHMM_model_path, n_perms, n_group_A, BIN_SIZE)
#
#
