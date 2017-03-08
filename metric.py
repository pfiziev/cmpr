from itertools import izip
import math
from utils import *

__author__ = 'Fiziev'

OVERALL_SCORE = 'overall_score'
BINS = 'bins'

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

def learn_metric_from_all_replicates(segmentations, states, BIN_SIZE):
    echo('Using all replicates to learn the metric')
    PSEUDO_COUNT = 1

    metrics = {}
    for dataset in segmentations:
        metric = dict((s1, dict((s2, PSEUDO_COUNT) for s2 in states)) for s1 in states)
        for rep_dataset in segmentations:
            if rep_dataset == dataset:
                continue

            seg1 = segmentations[dataset]
            seg2 = segmentations[rep_dataset]
            for chrom in seg1:
                for s1, s2 in izip(iterstate(seg1[chrom], BIN_SIZE), iterstate(seg2[chrom], BIN_SIZE)):
                    metric[s1][s2] += 1

        for s1 in states:
            total = sum(metric[s1].values())
            for s2 in states:
                metric[s1][s2] /= float(total)

        metrics[dataset] = metric

    return metrics

def get_prob_vector(bins, datasets, metric, states):

    probs = dict((s, 0) for s in states)

    for s in states:
        for dataset, d_state in izip(datasets, bins):
            probs[s] += metric[dataset][d_state][s]

    total = sum(probs[s] for s in states)
    for s in states:
        probs[s] /= total

    return probs


def get_prob_vector_based_on_posteriors(posteriors, datasets, metric):

    states = sorted(metric.values()[0].keys(), key=state_index)

    probs = dict((s, 0) for s in states)

    for dataset, (dataset_posteriors) in izip(datasets, posteriors):
        for state, state_posterior in izip(states, dataset_posteriors):
            for s2 in states:
                probs[s2] += state_posterior * metric[dataset][state][s2]

    total = sum(probs[s] for s in states)

    for s in states:
        probs[s] /= total

    return probs


def KL(P, Q):
    return sum(P[s] * math.log(P[s] / Q[s], 2) for s in P)


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
    for dataset in sorted(metric_X):
        print dataset
        metric = metric_X[dataset]
        states = sorted(metric, key=state_key)
        print '\t'.join(['State'] + states + ['Entropy'])
        for s1 in states:
            print '\t'.join([s1] + ['%.2lf' % (100 * metric[s1][s2]) for s2 in states] + ['%.2lf' % (entropy(metric[s1]))])


def print_average_metric(states, metric_A, metric_B):
    print '\nRaw'
    print '\t'.join(['State'] + states)
    for s1 in states:
        a_probs = get_prob_vector([s1 for _ in metric_A], metric_A.keys(), metric_A, states)
        print '\t'.join([s1] + [str(symmetric_KL_divergence(
                                    a_probs,
                                    get_prob_vector([s2 for _ in metric_B], metric_B.keys(), metric_B, states)))
                                for s2 in states])
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

def generate_pure_state_transitions(states, metric_A, metric_B):

    pure_state_transitions = {}
    for s1 in states:
        pure_state_transitions[s1] = {}
        a_probs = get_prob_vector([s1 for _ in metric_A], metric_A.keys(), metric_A, states)
        for s2 in states:
            b_probs = get_prob_vector([s2 for _ in metric_B], metric_B.keys(), metric_B, states)

            pure_state_transitions[s1][s2] = {
                                              # STATE_DISTRIBUTION: flat_dict({GROUP_A: dict((s,
                                              #                                               p if s == s1 else
                                              #                                               (1 - p) / (len(states) - 1))
                                              #                                              for s in states),
                                              #                                GROUP_B: dict((s,
                                              #                                               p if s == s2 else
                                              #                                               (1 - p) / (len(states) - 1))
                                              #                                              for s in states)}),
                                              STATE_DISTRIBUTION: flat_dict({GROUP_A: a_probs,
                                                                             GROUP_B: b_probs}),
                                              REAL_COUNTS: PSEUDO_COUNT,
                                              REAL_FREQUENCY: None,
                                              SHUFFLED_FREQUENCY: None,
                                              P_VALUE: None}
    return pure_state_transitions

def print_pure_state_transitions(pure_state_transitions):
    states = sorted(pure_state_transitions, key=state_key)
    for key in [REAL_COUNTS, SHUFFLED_FREQUENCY, REAL_FREQUENCY, P_VALUE]:
        print key
        print '\t'.join(['State'] + states)
        for s1 in sorted(pure_state_transitions, key=state_key):
            print '\t'.join([s1] + map(str, [pure_state_transitions[s1][s2][key] for s2 in states]))
        print


flat_dict = lambda d: dict((g + ' ' + k, d[g][k]) for g in d for k in d[g])
