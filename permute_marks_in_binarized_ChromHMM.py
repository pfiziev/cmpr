import argparse
import gzip
from itertools import izip
import os
import shutil
import sys, random
import tempfile
import operator
import re
from chrom_compare import get_overall_and_per_state_diff_score
from metric import *
from utils import echo, read_segmentations, filter_chroms, smooth
from subprocess import call

CHROM_HMM_EXEC = '/Users/Fiziev/software/ChromHMM/ChromHMMp.jar'


def read_binarized(fname, random_chunks):
    echo('Reading:', fname)
    with gzip.open(fname) if fname.endswith('.gz') else open(fname) as in_f:
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

    return marks, ct, chrom, c_data

def clean_up_directory(dirname):
    for f in os.listdir(dirname):
        os.remove(os.path.join(dirname, f))

def compute_background_scores_by_shuffling_marks(bin_directory,
                                                 seg_fnames,
                                                 ChromHMM_model_path,
                                                 random_chunks,
                                                 n_group_A,
                                                 BIN_SIZE,
                                                 states,
                                                 n_perms=100,
                                                 compute_per_state_scores=False,
                                                 max_min_window=0,
                                                 to_smooth=False):

    marks = None
    chrom_lengths = {}
    ct_chrom_fnames = {}
    celltypes = set(re.sub(r'_(\d+)$', '', os.path.split(sf)[1].split('_segments')[0]) for sf in seg_fnames)

    out_dir = tempfile.mkdtemp()
    echo('temp dir:', out_dir)

    data = {}
    marks = None
    background_model = dict((k, {}) for k in [OVERALL_SCORE] + (states if compute_per_state_scores else []))

    for fname in os.listdir(bin_directory):
        if fname.endswith('_binary.txt') and any(ct in fname for ct in celltypes):
            fname = os.path.join(bin_directory, fname)
            c_marks, ct, original_chrom, c_data = read_binarized(fname, random_chunks)

            if c_data is None:
                continue

            if marks is None:
                marks = c_marks

            if ct not in data:
                data[ct] = {}
                ct_chrom_fnames[ct] = {}

            for chrom in c_data:
                data[ct][chrom] = c_data[chrom]
                ct_chrom_fnames[ct][chrom] = fname.replace('_' + original_chrom + '_', '_' + chrom + '_')
                chrom_lengths[chrom] = len(c_data[chrom][marks[0]])

    mark_celltype_order = dict((m, list(data.keys())) for m in marks)

    seen = set()
    # seen.add(tuple([d for m in mark_celltype_order for d in mark_celltype_order[m]]))
    for perm_idx in xrange(n_perms):
        if perm_idx >= 0:
            key = None
            while key is None or key in seen:
                for m in mark_celltype_order:
                    random.shuffle(mark_celltype_order[m])
                key = tuple([d for m in mark_celltype_order for d in mark_celltype_order[m]])
            seen.add(key)
        echo('combo:', perm_idx)
        for chrom in chrom_lengths:
            for ct_idx, ct in enumerate(ct_chrom_fnames):
                out_fname = re.sub(r'\.gz$', '', os.path.join(out_dir, os.path.split(ct_chrom_fnames[ct][chrom])[1]))

                echo('marks:', marks)
                echo('ct order:', [mark_celltype_order[m][ct_idx] for m in marks])

                echo('Writing to:', out_fname)
                with open(out_fname, 'w') as out_f:

                    out_f.write('\t'.join([ct, chrom]) + '\n')
                    out_f.write('\t'.join(marks) + '\n')

                    for bin_idx in xrange(chrom_lengths[chrom]):
                        out_f.write(
                            '\t'.join([data[mark_celltype_order[m][ct_idx]][chrom][m][bin_idx] for m in marks]) + '\n')
        echo('Running ChromHMM')
        call(['java', '-mx8000M', '-jar', CHROM_HMM_EXEC, 'MakeSegmentation', ChromHMM_model_path, out_dir, out_dir])

        process_shuffled_segmentations(n_group_A,
                                       out_dir,
                                       background_model,
                                       BIN_SIZE,
                                       compute_per_state_scores,
                                       max_min_window,
                                       to_smooth)

        clean_up_directory(out_dir)

    echo('Purging temp directory:', out_dir)
    shutil.rmtree(out_dir)

    # convert the score counts to a CDF
    for score_type in background_model:
        total_regions = sum(background_model[score_type].itervalues())
        total_so_far = 0

        for s in sorted(background_model[score_type], reverse=True):
            total_so_far += background_model[score_type][s]
            background_model[score_type][s] = total_so_far / float(total_regions)

        background_model[score_type] = sorted(background_model[score_type].items())

    return dict((score_type, sorted(background_model[score_type].items())) for score_type in background_model)


def process_shuffled_segmentations(n_group_A,
                                   out_dir,
                                   background_model,
                                   BIN_SIZE,
                                   compute_per_state_scores,
                                   max_min_window,
                                   to_smooth=False):

    seg_fnames = sorted([os.path.join(out_dir, f) for f in os.listdir(out_dir) if f.endswith('_segments.bed')])

    group_A_seg_fnames = seg_fnames[:n_group_A]
    group_B_seg_fnames = seg_fnames[n_group_A:]

    group_A_segmentations, states = read_segmentations(group_A_seg_fnames)
    group_B_segmentations, _ = read_segmentations(group_B_seg_fnames)

    chromosomes = sorted(reduce(operator.and_, [set(s.keys())
                                             for segs in [group_A_segmentations, group_B_segmentations]
                                             for s in segs.values()]),
                         key=lambda c: int(c.replace('chr','')) if re.search(r'^chr(\d+)$', c) else 100)
    # print chromosomes
    group_A_segmentations = filter_chroms(group_A_segmentations, chromosomes)
    group_B_segmentations = filter_chroms(group_B_segmentations, chromosomes)

    metric_A = learn_metric_from_all_replicates(group_A_segmentations, states, BIN_SIZE)
    # print_metric(metric_A)

    metric_B = learn_metric_from_all_replicates(group_B_segmentations, states, BIN_SIZE)
    # print_metric(metric_B)
    # print_average_metric(states, metric_A, metric_B)

    for chrom in chromosomes:
        echo('Chromosome:', chrom)

        chrom_segmentations_A = dict((d, chrom_segmentation_to_list(group_A_segmentations[d][chrom], BIN_SIZE)) for d in group_A_seg_fnames)
        chrom_segmentations_B = dict((d, chrom_segmentation_to_list(group_B_segmentations[d][chrom], BIN_SIZE)) for d in group_B_seg_fnames)


        chunk_scores = get_overall_and_per_state_diff_score(chrom,
                                                            chrom_segmentations_A,
                                                            metric_A,
                                                            chrom_segmentations_B,
                                                            metric_B,
                                                            BIN_SIZE,
                                                            states,
                                                            None,
                                                            compute_per_state_scores,
                                                            max_min_window,
                                                            background_chunk=True)
        if to_smooth:
            chunk_scores = dict((score_type, smooth(chunk_scores[score_type])) for score_type in chunk_scores)

        for score_type in background_model:
            for s in chunk_scores[score_type]:
                abs_score = abs(s)
                if abs_score not in background_model[score_type]:
                    background_model[score_type][abs_score] = 0

                background_model[score_type][abs_score] += 1





if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', nargs='+', help='directory with binarized chromosome files')
    parser.add_argument('-s', nargs='+', help='original chromhmm segmentations')

    parser.add_argument('-n', help='n permutations', type=int)
    parser.add_argument('-a', help='number of samples in group A', type=int)
    parser.add_argument('--bin-size', default=200, help='bin size', type=int)
    parser.add_argument('-m', help='ChromHMM model')

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        exit(1)

    n_perms = args.n
    bin_directory = args.i
    seg_fnames = args.s
    ChromHMM_model_path = args.m
    n_group_A = args.a
    BIN_SIZE = args.bin_size

    # compute_background_scores_by_shuffling_marks(bin_directory, seg_fnames, ChromHMM_model_path, n_perms, n_group_A, BIN_SIZE)


