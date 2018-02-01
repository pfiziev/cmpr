from itertools import izip
import math
from utils import *

__author__ = 'Fiziev'

OVERALL_SCORE = 'overall_score'
BINS = 'bins'
CONSTITUTIVE = 'constitutive'

def _learn_metric(segmentations, states, BIN_SIZE, other_group=None):
    PSEUDO_COUNT = 1

    metrics = {}
    for dataset in segmentations:
        metric = dict((s1, dict((s2, PSEUDO_COUNT) for s2 in states)) for s1 in states)
        closest_distance, closest = find_closest(segmentations, dataset, segmentations[dataset], states, BIN_SIZE)
        print 'Dataset:', dataset, 'Closest:', closest, 'Distance:', closest_distance
        if other_group is not None:
            closest_distance_from_other_group, closest_from_other_group = find_closest(other_group,
                                                                                       dataset,
                                                                                       segmentations[dataset],
                                                                                       states,
                                                                                       BIN_SIZE)
            if closest_distance > closest_distance_from_other_group:
                echo('WARNING: Dataset has a closer neighbor in the other group!')
                echo('WARNING: Dataset:',
                     dataset,
                     'closest in other group:',
                     closest_from_other_group,
                     'Distance:',
                     closest_distance_from_other_group)
            else:
                echo('No closer datasets in the other group')

        seg1 = segmentations[dataset]
        seg2 = segmentations[closest]
        for chrom in seg1:
            for s1, s2 in izip(iterstate(seg1[chrom], BIN_SIZE), iterstate(seg2[chrom], BIN_SIZE)):
                metric[s1][s2] += 1

        for s1 in states:
            total = sum(metric[s1].values())
            for s2 in states:
                metric[s1][s2] /= float(total)

        metrics[dataset] = metric

    return metrics


metric_cache = {}
def clear_metric_cache():
    for k in list(metric_cache):
        del metric_cache[k]


def learn_metric_from_all_replicates(segmentations, states, BIN_SIZE, posteriors_dir=None, clear_cache=False):

    if clear_cache:
        clear_metric_cache()

    if posteriors_dir is not None:
        return learn_metric_from_all_replicates_with_posteriors(segmentations, states, posteriors_dir)

    # echo('Using all replicates to learn the metric')

    # The pseudo count assumes that each state pair occurred at least once when comparing each pair of datasets
    PSEUDO_COUNT = len(segmentations) - 1

    metrics = {}
    for dataset in segmentations:
        metric = dict((s1, dict((s2, PSEUDO_COUNT) for s2 in states)) for s1 in states)
        for rep_dataset in segmentations:
            if rep_dataset == dataset:
                continue

            key = (dataset, rep_dataset)
            if key not in metric_cache:
                _metric = dict((s1, dict((s2, 0) for s2 in states)) for s1 in states)

                seg1 = segmentations[dataset]
                seg2 = segmentations[rep_dataset]
                for chrom in seg1:
                    for s1, s2 in izip(iterstate(seg1[chrom], BIN_SIZE), iterstate(seg2[chrom], BIN_SIZE)):
                        _metric[s1][s2] += 1
                metric_cache[key] = _metric

            for s1 in metric:
                for s2 in metric:
                    metric[s1][s2] += metric_cache[key][s1][s2]

        for s1 in states:
            total = sum(metric[s1].values())
            for s2 in states:
                metric[s1][s2] /= float(total)

        metrics[dataset] = metric

    return metrics


def compute_average_metric(metric_A, metric_B, states):
    average_metric = dict((s1, dict((s2, 0) for s2 in states)) for s1 in states)

    total_datasets = len(metric_A) + len(metric_B)

    for m in [metric_A, metric_B]:
        for d in m:
            for s1 in m[d]:
                for s2 in m[d][s1]:
                    average_metric[s1][s2] += m[d][s1][s2]

    average_metric = dict((s1, dict((s2, average_metric[s1][s2] / total_datasets) for s2 in states)) for s1 in states)
    return dict((d, average_metric) for d in metric_A), dict((d, average_metric) for d in metric_B)


#
# def learn_metric_from_all_replicates_test(segmentations, states, BIN_SIZE, posteriors_dir=None, clear_cache=False):
#
#     metrics = learn_metric_from_all_replicates(segmentations, states, BIN_SIZE)
#
#     def _slice(seg, chunk_start, chunk_end):
#         sliced_seg = []
#         i = 0
#         while not (seg[i][0] <= chunk_start < seg[i][1]):
#             i += 1
#         offset = chunk_start
#         while i < len(seg) and overlap(chunk_start, chunk_end, seg[i][0], seg[i][1]) > 0:
#             s = max(chunk_start, seg[i][0])
#             e = min(chunk_end, seg[i][1])
#             sliced_seg.append([s - offset, e - offset, seg[i][2]])
#             i += 1
#
#         return sliced_seg
#
#
#     for n_chunks in xrange(0, 200, 20):
#         if n_chunks == 0:
#             n_chunks = 1
#         all_diffs = []
#         for i in xrange(100):
#             chromosomes = sorted(segmentations.values()[0])
#             random_chunks = generate_random_chunks(segmentations, chromosomes, BIN_SIZE, N_RANDOM_CHUNKS=n_chunks)
#             chunk_segs = dict((d, dict((chrom + '_' + str(chunk_no),
#                                         _slice(segmentations[d][chrom], chunk_start * BIN_SIZE, chunk_end * BIN_SIZE))
#                                             for chrom in random_chunks
#                                                 for chunk_no, (chunk_start, chunk_end) in enumerate(random_chunks[chrom])))
#                               for d in segmentations)
#
#             c_metrics = learn_metric_from_all_replicates(chunk_segs, states, BIN_SIZE, clear_cache=True)
#             diff = 0
#             for d in metrics:
#                 full_metric = metrics[d]
#                 partial_metric = c_metrics[d]
#                 for s1 in states:
#                     for s2 in states:
#                         if abs(full_metric[s1][s2] - partial_metric[s1][s2]) >= 0.1:
#                             diff += 1
#             diff /= float(len(metrics))
#             all_diffs.append(diff)
#
#         print 'n_chunks=', n_chunks, mean_and_std(all_diffs), min(all_diffs), max(all_diffs)
#
#     return metrics



def learn_metric_from_all_replicates_with_posteriors(segmentations, states, posteriors_dir):
    echo('Using all replicates to learn the metric with posteriors')
    PSEUDO_COUNT = len(segmentations) - 1

    metrics = {}
    for dataset in segmentations:
        echo(dataset)
        metric = dict((s1, dict((s2, PSEUDO_COUNT) for s2 in states)) for s1 in states)
        for rep_dataset in segmentations:
            if rep_dataset == dataset:
                continue

            seg1 = segmentations[dataset]

            for chrom in seg1:
                echo(chrom)
                seg1_posteriors = iter_posteriors(get_chrom_posteriors_fname(posteriors_dir, dataset, chrom))
                seg2_posteriors = iter_posteriors(get_chrom_posteriors_fname(posteriors_dir, rep_dataset, chrom))

                for s1_probs, s2_probs in izip(seg1_posteriors, seg2_posteriors):
                    s1_state_idx = filter(lambda s_idx: s1_probs[s_idx] >= MIN_POSTERIOR_TO_CONSIDER, range(len(states)))
                    s2_state_idx = filter(lambda s_idx: s2_probs[s_idx] >= MIN_POSTERIOR_TO_CONSIDER, range(len(states)))

                    for s1_i in s1_state_idx:
                        for s2_i in s2_state_idx:
                            metric[states[s1_i]][states[s2_i]] += s1_probs[s1_i] * s2_probs[s2_i]

                    # for s1_i, s1 in enumerate(states):
                    #     for s2_i, s2 in enumerate(states):
                    #         metric[s1][s2] += s1_probs[s1_i] * s2_probs[s2_i]

        for s1 in states:
            total = sum(metric[s1].values())
            for s2 in states:
                metric[s1][s2] /= float(total)

        metrics[dataset] = metric

    return metrics


def iter_posteriors(fname):
    with open_file(fname) as f:
        f.readline()
        f.readline()
        for line in f:
            yield tuple(map(lambda v: round(float(v), 3), line.strip().split()))


def read_posteriors_as_list(fname):
    return list(iter_posteriors(fname))

def get_celltype(fname):
    return re.sub(r'_(\d+)(_coreMarks)?$', '', os.path.split(fname)[1].split('_segments')[0])


def get_chrom_posteriors_fname(posteriors_dir, fname, chrom):

    prefix = os.path.split(fname)[1].split('_segments')[0]
    uncompressed_fname = os.path.join(posteriors_dir, prefix + '_' + chrom + '_posterior.txt')
    if os.path.isfile(uncompressed_fname):
        return uncompressed_fname
    else:
        return uncompressed_fname + '.gz'


def get_prob_vector(bins, datasets, metric, states):

    probs = dict((s, 0) for s in states)

    for s in states:
        for dataset, d_state in izip(datasets, bins):
            probs[s] += metric[dataset][d_state][s]

    total = sum(probs[s] for s in states)
    for s in states:
        probs[s] /= total

    return probs

MIN_POSTERIOR_TO_CONSIDER = 0.05

def get_prob_vector_based_on_posteriors(posteriors, datasets, metric, states):

    probs = dict((s, 0) for s in states)

    for dataset, (dataset_posteriors) in izip(datasets, posteriors):
        for state, state_posterior in izip(states, dataset_posteriors):
            if state_posterior >= MIN_POSTERIOR_TO_CONSIDER :
                for s2 in states:
                    probs[s2] += state_posterior * metric[dataset][state][s2]

    total = sum(probs[s] for s in states)

    for s in states:
        probs[s] /= total

    return probs


def KL(P, Q):
    kl = sum(P[s] * math.log(P[s] / Q[s], 2) for s in P)
    return kl if kl > 0 else 0  # protect from underflow errors that can yield a negative number


def symmetric_KL_divergence(probs_a, probs_b):
    return (KL(probs_a, probs_b) + KL(probs_b, probs_a)) / 2.


# def hellinger_distance(probs_a, probs_b):
#     return math.sqrt(1 - sum(math.sqrt(probs_a[s] * probs_b[s]) for s in probs_a))

def hellinger_distance(probs_a, probs_b):
    states = sorted(probs_a)
    return eucld([math.sqrt(probs_a[s]) for s in states],
                 [math.sqrt(probs_b[s]) for s in states]) / math.sqrt(2)

# KL = hellinger_distance
# KL = lambda P, Q: max(abs(P[x] - Q[x]) for x in P)

def state_index(s):
    return int(s[1:]) if not s.endswith('NONE') else 0


def find_closest(segmentations, dataset, dataset_segmentation, states, BIN_SIZE):
    return min([(log_enrichment_distance(segmentations[d],
                                         dataset_segmentation,
                                         states,
                                         BIN_SIZE),
                 d) for d in segmentations if d != dataset])

    # return min([d for d in segmentations if d != dataset],
    #            key=lambda _d: log_enrichment_distance(segmentations[_d],
    #                                                   segmentations[dataset],
    #                                                   states,
    #                                                   BIN_SIZE))


def iterstate(chrom_segmentation, BIN_SIZE):
    for start, end, state in chrom_segmentation:
        for bin_start in xrange(start, end, BIN_SIZE):
            yield state


def chrom_segmentation_to_list(chrom_segmentation, BIN_SIZE):
    return list(iterstate(chrom_segmentation, BIN_SIZE))


def log_enrichment_distance(seg1, seg2, states, BIN_SIZE):

    PSEUDO_COUNT = 1

    overlap = dict((s1, dict((s2, 0) for s2 in states)) for s1 in states)

    # compute the overlap between every pair of states
    for chrom in seg1:
        if chrom in seg2:
            for s1, s2 in izip(iterstate(seg1[chrom], BIN_SIZE), iterstate(seg2[chrom], BIN_SIZE)):
                overlap[s1][s2] += 1

    # compute state lengths in cell type 1
    state_length_1 = dict((s, sum(overlap[s].values())) for s in states)

    # compute state lengths in cell type 2
    state_length_2 = dict((s, sum(overlap[s1][s] for s1 in states)) for s in states)

    # compute the genome length
    genome_length = sum(state_length_1.values())

    # compute enrichments for the overlap of each state with itself and add a pseudo-count of 1
    diag_enrichments = [float(PSEUDO_COUNT + overlap[s][s]) /
                        (PSEUDO_COUNT + float(state_length_1[s] * state_length_2[s]) / genome_length)
                        for s in states]

    # print [math.log(e, 2) for e in diag_enrichments]
    # print [overlap[s][s] for s in states]
    # return the negative log of the geometric mean
    return -sum([math.log(e, 2) for e in diag_enrichments])

def entropy(d):
    e = 0
    for p in d.values():
        if p > 0:
            e -= p * math.log(p, 2)
    return e

def print_metric(metric_X):
    to_print = ''
    for dataset in sorted(metric_X):
        to_print += dataset + '\n'
        metric = metric_X[dataset]
        states = sorted(metric, key=state_key)
        to_print += '\t'.join(['State'] + states + ['Entropy']) + '\n'
        for s1 in states:
            to_print += '\t'.join([s1] + ['%.2lf' % (100 * metric[s1][s2]) for s2 in states] + ['%.2lf' % (entropy(metric[s1]))]) + '\n'

    echo(to_print)


def print_consistent_state_transition_scores(states, metric_A, metric_B):
    to_print = '\nRaw' + '\n'
    to_print += '\t'.join(['State'] + states) + '\n'
    for s1 in states:
        a_probs = get_prob_vector([s1 for _ in metric_A], metric_A.keys(), metric_A, states)
        to_print += '\t'.join([s1] + [str(symmetric_KL_divergence(
                                    a_probs,
                                    get_prob_vector([s2 for _ in metric_B], metric_B.keys(), metric_B, states)))
                                for s2 in states]) + '\n'

    echo(to_print)

    # print '\nSquares'
    # print '\t'.join(['State'] + states)
    # for s1 in states:
    #     a_probs = get_prob_vector([(0, s1) for _ in metric_A], metric_A.keys(), metric_A, states)
    #     print '\t'.join([s1] + [str(symmetric_KL_divergence(
    #                                 a_probs,
    #                                 get_prob_vector([(0, s2) for _ in metric_B], metric_B.keys(), metric_B, states)) ** 2)
    #                             for s2 in states])

GROUP_A = 'A'
GROUP_B = 'B'
STATE_DISTRIBUTION = 'state distribution'

REAL_COUNTS = 'real counts'
REAL_FREQUENCY = 'real frequency'

SHUFFLED_FREQUENCY = 'shuffled frequency'
SHUFFLED_COUNTS = 'shuffled counts'

P_VALUE = 'p-value'
PSEUDO_COUNT = 10

# def generate_pure_state_transitions(states, metric_A, metric_B):
#
#
#     pure_state_transitions = {}
#     for s1 in states:
#         pure_state_transitions[s1] = {}
#         a_probs = get_prob_vector([s1 for _ in metric_A], metric_A.keys(), metric_A, states)
#         for s2 in states:
#             b_probs = get_prob_vector([s2 for _ in metric_B], metric_B.keys(), metric_B, states)
#
#             pure_state_transitions[s1][s2] = {
#                                               # STATE_DISTRIBUTION: flat_dict({GROUP_A: dict((s,
#                                               #                                               p if s == s1 else
#                                               #                                               (1 - p) / (len(states) - 1))
#                                               #                                              for s in states),
#                                               #                                GROUP_B: dict((s,
#                                               #                                               p if s == s2 else
#                                               #                                               (1 - p) / (len(states) - 1))
#                                               #                                              for s in states)}),
#                                               STATE_DISTRIBUTION: flat_dict({GROUP_A: a_probs,
#                                                                              GROUP_B: b_probs}),
#                                               REAL_COUNTS: PSEUDO_COUNT,
#                                               REAL_FREQUENCY: None,
#                                               SHUFFLED_FREQUENCY: None,
#                                               P_VALUE: None}
#     return pure_state_transitions


def generate_consistent_state_probabilities(states, metric_A, metric_B):

    consistent_state_probs = {GROUP_A: {},
                              GROUP_B: {}}

    for grp_id, metric in [(GROUP_A, metric_A), (GROUP_B, metric_B)]:
        for s in states:
            consistent_state_probs[grp_id][s] = get_prob_vector([s for _ in metric], metric.keys(), metric, states)

    return consistent_state_probs

#
# def print_pure_state_transitions(pure_state_transitions):
#     states = sorted(pure_state_transitions, key=state_key)
#     for key in [REAL_COUNTS, SHUFFLED_FREQUENCY, REAL_FREQUENCY, P_VALUE]:
#         print key
#         print '\t'.join(['State'] + states)
#         for s1 in sorted(pure_state_transitions, key=state_key):
#             print '\t'.join([s1] + map(str, [pure_state_transitions[s1][s2][key] for s2 in states]))
#         print


flat_dict = lambda d: dict((g + ' ' + k, d[g][k]) for g in d for k in d[g])


def new_score_dict(states, compute_per_state_scores):
    return dict((k, {}) for k in [OVERALL_SCORE] + (states if compute_per_state_scores else []))
