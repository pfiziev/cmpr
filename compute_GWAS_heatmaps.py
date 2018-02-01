import argparse
import os, sys
from utils import *


def filter_files(files):
    # return [f for f in files if 'smooth' in f]
    return [f for f in files if 'raw' in f] if args.raw else [f for f in files if 'smooth' in f]

parser = argparse.ArgumentParser()

parser.add_argument('-d', help='directory with differential scores')
parser.add_argument('-c', help='directory with hard state counts')
parser.add_argument('-s', help='directory with soft counts')
parser.add_argument('-a', help='directory with one vs all scores')
parser.add_argument('-t', help='directory with constitutive scores')

parser.add_argument('--raw', action='store_true', default=False, help='use the raw scores, by default use smoothed')

parser.add_argument('-n', type=int, help='number of states')
parser.add_argument('-o', '--output', help='output prefix')

args = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    exit(1)


HARD_COUNTS = 'hard_counts'
CONSTITUTIVE = 'constitutive'

DIFF_SCORE_BEST = 'diff_score_best'
DIFF_SCORE_MEDIAN = 'diff_score_median'

SOFT_COUNTS = 'soft_counts'
ONE_VS_ALL_SCORE = 'one_vs_all_score'

OVERALL_SCORE = 'overall_score'
FIRST_CELLTYPE = 1
SECOND_CELLTYPE = 2

score_types = [HARD_COUNTS, DIFF_SCORE_BEST, DIFF_SCORE_MEDIAN, SOFT_COUNTS, ONE_VS_ALL_SCORE, CONSTITUTIVE]
# score_types = [HARD_COUNTS, DIFF_SCORE_BEST, DIFF_SCORE_MEDIAN, SOFT_COUNTS]

states = ['E' + str(s) for s in range(1, args.n + 1)]

# data = dict((st, dict((s, {}) for s in states + ([OVERALL_SCORE] if st in [DIFF_SCORE_BEST, DIFF_SCORE_MEDIAN] else [])))
#             for st in score_types)

filenames = dict((score_type, dict((s, {}) for s in states))
                   for score_type in score_types)

echo('Processing files in:', args.d)
celltypes = set()

for fname in filter_files(os.listdir(args.d)):
    full_fname = os.path.join(args.d, fname)

    parts = fname.replace('gwas_', '').split('.')
    ct1, ct2 = parts[0].replace('-', '_'), parts[2].replace('-', '_')
    celltypes.add(ct1)
    celltypes.add(ct2)
    state = parts[-3]
    if state in states:
        for score_type in [DIFF_SCORE_MEDIAN, DIFF_SCORE_BEST]:

            if ct1 not in filenames[score_type][state]:
                filenames[score_type][state][ct1] = []

            if ct2 not in filenames[score_type][state]:
                filenames[score_type][state][ct2] = []

            filenames[score_type][state][ct1].append((full_fname, FIRST_CELLTYPE))
            filenames[score_type][state][ct2].append((full_fname, SECOND_CELLTYPE))
    elif state == OVERALL_SCORE:
        pass

echo('Processing files in:', args.s)
for fname in filter_files(os.listdir(args.s)):
    full_fname = os.path.join(args.s, fname)

    parts = fname.replace('gwas_', '').split('.')
    celltype = parts[0].replace('-', '_')
    celltypes.add(celltype)
    state = parts[-3]
    if state in states:
        filenames[SOFT_COUNTS][state][celltype] = full_fname

echo('Processing files in:', args.c)
for fname in os.listdir(args.c):
    full_fname = os.path.join(args.c, fname)

    parts = fname.replace('gwas_', '').split('.')
    ct = parts[0].replace('-', '_')
    celltypes.add(ct)

    state = parts[2]
    if state in states:
        filenames[HARD_COUNTS][state][ct] = full_fname

echo('Processing files in:', args.a)
for fname in filter_files(os.listdir(args.a)):
    full_fname = os.path.join(args.a, fname)

    parts = fname.replace('gwas_', '').split('.')
    ct = parts[2].replace('-', '_')
    celltypes.add(ct)
    state = parts[-3]
    if state in states:
        filenames[ONE_VS_ALL_SCORE][state][ct] = full_fname

echo('Processing files in:', args.t)
for fname in filter_files(os.listdir(args.t)):
    full_fname = os.path.join(args.t, fname)

    parts = fname.replace('gwas_', '').split('.')

    state = parts[-3]
    if state in states:
        filenames[CONSTITUTIVE][state] = full_fname

echo('Celltypes:', sorted(celltypes))


def read_gwas_stats(fname, celltype_number):
    echo('Reading:', fname)

    stats = {}
    with gzip.open(fname) as in_f:
        for l in in_f:
            buf = l.strip().split('\t')
            pval = float(buf[0])

            if pval == 0:
                pval = 10 ** -20

            pubmed = buf[3]
            trait = buf[4]
            avg_bgr_idx = float(buf[5])
            avg_fgr_idx = float(buf[6])

            key = trait + ', p:' + pubmed

            if (celltype_number == FIRST_CELLTYPE and avg_fgr_idx < avg_bgr_idx) or \
               (celltype_number == SECOND_CELLTYPE and avg_fgr_idx > avg_bgr_idx):
                sign = 1
            else:
                sign = -1

            stats[key] = sign * pval

            # if pval <= P_VAL_THRESHOLD:
            #     if (celltype_number == FIRST_CELLTYPE and avg_fgr_idx < avg_bgr_idx) or \
            #        (celltype_number == SECOND_CELLTYPE and avg_fgr_idx > avg_bgr_idx):
            #         key = trait + ', p:' + pubmed
            #         stats[key] = pval
            #

    return stats


all_traits = {}
all_heatmaps = {}

for score_type in score_types:
    echo(score_type)
    all_heatmaps[score_type] = {}
# for score_type in [HARD_COUNTS]:
    for state in sorted(filenames[score_type]):
        heatmap = {}

        if state not in all_traits:
            all_traits[state] = {}

        if score_type == CONSTITUTIVE:
            all_heatmaps[score_type][state] = {CONSTITUTIVE:
                                                   read_gwas_stats(filenames[score_type][state], FIRST_CELLTYPE)}

            for trait in all_heatmaps[score_type][state][CONSTITUTIVE]:
                    if trait not in all_traits[state]:
                        all_traits[state][trait] = 1

                    all_traits[state][trait] = min(all_traits[state][trait],
                                                   abs(all_heatmaps[score_type][state][CONSTITUTIVE][trait]))

            continue

        celltypes = sorted(filenames[score_type][state])
        for celltype in celltypes:

            if score_type in [HARD_COUNTS, SOFT_COUNTS, ONE_VS_ALL_SCORE]:
                heatmap[celltype] = read_gwas_stats(filenames[score_type][state][celltype], FIRST_CELLTYPE)


            elif score_type in [DIFF_SCORE_BEST, DIFF_SCORE_MEDIAN]:
                all_celltype_state_states = {}
                for fname, celltype_number in filenames[score_type][state][celltype]:
                    stats = read_gwas_stats(fname, celltype_number)
                    for key in stats:
                        if key not in all_celltype_state_states:
                            all_celltype_state_states[key] = []
                        all_celltype_state_states[key].append(stats[key])

                # use min for CELLTYPE_SCORES and DIFF_SCORE_BEST, else median

                def _min(array):
                    a = [(abs(p), sign(p)) for p in array]
                    min_v, min_s = min(a)
                    return min_v * min_s

                def _median(array):
                    a = [(abs(p), sign(p)) for p in array]
                    min_v, min_s = sorted(a)[len(a) / 2]
                    return min_v * min_s

                get_value = _min if score_type in [DIFF_SCORE_BEST] else _median

                heatmap[celltype] = dict((key, get_value(all_celltype_state_states[key])) for key in all_celltype_state_states)

            for trait in heatmap[celltype]:
                    if trait not in all_traits[state]:
                        all_traits[state][trait] = 1

                    all_traits[state][trait] = min(all_traits[state][trait], abs(heatmap[celltype][trait]))

        all_heatmaps[score_type][state] = heatmap

P_VAL_THRESHOLD = 10 ** -3.5
all_traits = dict((state,
                   sorted(t for t in all_traits[state] if all_traits[state][t] <= P_VAL_THRESHOLD))
                  for state in states)


for score_type in score_types:
    echo('outputting:', score_type)
    best_total = 0
    best_combo = None

    # for score_type in [HARD_COUNTS]:
    for state in sorted(all_heatmaps[score_type]):
        total = 0
        heatmap = all_heatmaps[score_type][state]
        out_fname = args.output + '.' + score_type + '.' + state + '.heatmap.txt'

        with open(out_fname, 'w') as out_f:
            _celltypes = [CONSTITUTIVE] if score_type == CONSTITUTIVE else celltypes

            out_f.write(out_fname.replace('.heatmap.txt', '') + ('\t' * len(_celltypes)) +'\n')
            out_f.write('\t'.join([''] + _celltypes) + '\n')

            for trait in sorted(all_traits[state]):
                p_values = []

                for celltype in _celltypes:
                    p_value = heatmap.get(celltype, {}).get(trait, 1)
                    p_values.append(-sign(p_value) * math.log(abs(p_value), 10))

                out_f.write('\t'.join([trait] + map(str, p_values)) + '\n')
                total += sum(p for p in p_values if p > 0)

        if total >= best_total:
            best_total = total
            best_combo = state

    echo('Best state for', score_type, ':', best_combo)

echo('Done')

