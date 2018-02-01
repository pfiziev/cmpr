"""
    Here is the place to put utility methods that are shared by the modules.
"""
import gzip
from itertools import izip
import json
import os
import re

import sys
import string
import datetime
import math
import ctypes
import random
import gc

random.seed(42)

ROOT_DIR = os.path.split(__file__)[0]


def error(msg):
    print >> sys.stderr, 'ERROR: %s' % msg
    exit(1)



global_stime = datetime.datetime.now()


def elapsed(message = None):
    """ Measures how much time has elapsed since the last call of this method and the beginning of the execution.
        If 'message' is given, the message is printed out with the times.
    """
    print "[Last: " + str(datetime.datetime.now() - elapsed.stime) + ', Elapsed time: '+str(datetime.datetime.now() - global_stime)+ "] %s" % message if message is not None else ""
    elapsed.stime = datetime.datetime.now()
elapsed.stime = datetime.datetime.now()


def open_log(fname):
    """ Opens a log file
    """
    open_log.logfile = open(fname, 'w', 1)
    open_log.start_time = datetime.datetime.now()

# def logm(message):
#     """ Logs a message with a time stamp to the log file opened by open_log
#     """
#     print "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), message)
#     open_log.logfile.write("[ %s ] %s\n" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), message))


def echo(*message):
    to_print = "[pid:%d %s] %s" % (os.getpid(),
                                 datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                 ' '.join(map(str, message)))
    print to_print
    if hasattr(open_log, 'logfile'):
        open_log.logfile.write(to_print + '\n')


def echo_to_stdder(*message):
    to_print = "[%s] %s" % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), ' '.join(map(str, message)))
    print >>sys.stderr, to_print
    if hasattr(open_log, 'logfile'):
        open_log.logfile.write(to_print + '\n')

def close_log():
    """ Closes the log file opened by open_log
    """
    echo('Elapsed time:', datetime.datetime.now() - open_log.start_time)
    open_log.logfile.close()

def mean(array):
    return float(sum(array)) / len(array) if len(array) > 0 else 0

def std(array):
    m = mean(array)
    return math.sqrt(sum((x - m) ** 2 for x in array) / float(len(array)))

def unbiased_std(array):
    m = mean(array)
    return math.sqrt(sum((x - m) ** 2 for x in array) / float(len(array) - 1))


def variance(array):
    return std(array) ** 2


def mean_and_std(array):
    m = mean(array)
    return m, math.sqrt(sum((x - m) ** 2 for x in array) / float(len(array)))

def mean_and_std_and_variance(array):
    m = mean(array)
    var = sum((x - m) ** 2 for x in array) / float(len(array))
    return m, math.sqrt(var), var

def standardize(array):
    m, s = mean_and_std(array)
    if s > 0:
        return [(v - m) / s for v in array]
    else:
        return [(v - m) for v in array]

def eucld(x, y):
    return math.sqrt(sum((xx - yy) ** 2 for xx, yy in izip(x, y)))

rmse = eucld


def R2(responses, predictions):
    mean_response = mean(responses)
    var_response = sum((r - mean_response) ** 2 for r in responses)

    return 1 - (float(sum((r - p) ** 2 for r, p in izip(responses, predictions))) / var_response)


def pearsonr(responses, predictions):
    mean_response = mean(responses)
    var_response = sum((r - mean_response) ** 2 for r in responses)

    mean_prediction = mean(predictions)
    var_prediction = sum((r - mean_prediction) ** 2 for r in predictions)

    return sum((r - mean_response) * (p - mean_prediction)
                   for r, p in izip(responses, predictions)) / math.sqrt(var_response * var_prediction)


def open_file1(fname, mode='r'):
    return gzip.open(fname, mode) if fname.endswith('.gz') else open(fname, mode)


import subprocess, io


class _GZipFileWriter:
    def __init__(self, fileName):
        # echo('Constructing _GZipFileWriter:', fileName)

        self.f = open(fileName, 'w')
        self.p = subprocess.Popen(['gzip'], stdin=subprocess.PIPE, stdout=self.f)
        # echo('Starting gzip child process:', self.p.pid)

    def __enter__(self):
        return self.p.stdin

    def __exit__(self, exc_type, exc_value, traceback):
        self.p.stdin.close()
        # echo('Waiting on', self.p.pid)
        self.p.wait()
        self.f.close()


def open_file(fname, mode='r'):
    if fname.endswith('.gz'):

        if mode == 'w':
            return _GZipFileWriter(fname)
        elif mode == 'r':
            return io.BufferedReader(gzip.open(fname, 'r'))
            # in_f = open(fname, 'rb')
            # p = subprocess.Popen(['gunzip'], stdout=subprocess.PIPE, stdin=in_f)
            # return p.stdout
        else:
            error('Only w and r are supported modes for gzipped files. fname=' + fname + ', mode=' + mode)

    else:
        return open(fname, mode)

#
# def reap_gzip_zombies():
#     to_remove = []
#     for i in xrange(len(reap_gzip_zombies.procs)):
#         p = reap_gzip_zombies.procs[i]
#         return_code = p.poll()
#
#         if return_code is not None:
#             echo('Reaping:', p, p.pid, return_code)
#             p.wait()
#             to_remove.append(i)
#
#     for i in to_remove:
#         reap_gzip_zombies.procs.pop(i)
#
# reap_gzip_zombies.procs = []

def append_and_unlink(from_fname, to_fname):
    # append the file content of from_file to to_file and delete from_file
    if type(to_fname) is str:
        to_f = open_file(to_fname, 'a')
    else:
        to_f = to_fname

    with open_file(from_fname) as from_f:
        to_f.write(from_f.read())

    os.unlink(from_fname)

    if type(to_fname) is str:
        to_f.close()


def not_None(thing):
    return thing is not None

def sign(value):
    return 0 if value == 0 else -1 if value < 0 else 1

def parse(s, types, splitter='\t'):
    buf = s.split(splitter)
    if len(buf) != len(types):
        error('Cannot parse string: ' + s + ' with types: ' + str(types))
    return [t(b) for t, b in zip(types, buf)]


def median(array):
    so = sorted(array)
    n = len(so)
    if n % 2:
        return so[n / 2]
    else:
        return float(so[n / 2 - 1] + so[n / 2]) / 2

def n_combinations(N, k):
    if (k > N) or (N < 0) or (k < 0):
        return 0L
    N, k = map(long,(N,k))
    top = N
    val = 1L
    while top > (N-k):
        val *= top
        top -= 1
    n = 1L
    while n < k+1L:
        val /= n
        n += 1
    return val

def box(data):
    return sum(data) / len(data)

def gaussian(data):
    mean = len(data) / 2
    std = float(len(data)) / 6
    if not hasattr(gaussian, 'n'):
        gaussian.n = sum(math.exp(-((i - mean)**2) / (2 * std ** 2)) for i, v in enumerate(data))
    return sum(math.exp(-((i - mean)**2) / (2 * std ** 2)) * v for i, v in enumerate(data)) / gaussian.n

PRECISION = 5

def smooth(wig_data, filter='gaussian', bandwidth=10):
    result = [0] * len(wig_data)
    window = [0] * (2 * bandwidth + 1)

    def copy_window(window, wig_data, bin_idx, bandwidth):
        for w_idx, score_idx in enumerate(xrange(bin_idx - bandwidth, bin_idx + bandwidth + 1)):
            if score_idx < 0 or score_idx >= len(wig_data):
                value = 0
            else:
                value = wig_data[score_idx]
            window[w_idx] = value

    for bin_idx in xrange(len(wig_data)):
        copy_window(window, wig_data, bin_idx, bandwidth)
        # window = wig_data[max(bin_idx - bandwidth, 0):
        #                   min(bin_idx + bandwidth + 2, len(wig_data))]
        if filter == 'box':
            result[bin_idx] = box(window)
        elif filter == 'gaussian':
            result[bin_idx] = gaussian(window)
        else:
            print >>sys.stderr, 'ERROR: unsupported filter:', filter
            exit(1)

        result[bin_idx] = round(result[bin_idx], PRECISION)

    return result


def smooth_in_place(wig_data, filter='gaussian', bandwidth=10):

    window = [0] * bandwidth + wig_data[:bandwidth + 1]

    def copy_window(window, wig_data, bin_idx, bandwidth):
        for i in xrange(len(window) - 1):
            window[i] = window[i + 1]
        window[-1] = wig_data[bin_idx + bandwidth] if bin_idx + bandwidth < len(wig_data) else 0

    for bin_idx in xrange(len(wig_data)):
        copy_window(window, wig_data, bin_idx, bandwidth)
        # window = wig_data[max(bin_idx - bandwidth, 0):
        #                   min(bin_idx + bandwidth + 2, len(wig_data))]
        if filter == 'box':
            wig_data[bin_idx] = box(window)
        elif filter == 'gaussian':
            wig_data[bin_idx] = gaussian(window)
        else:
            print >>sys.stderr, 'ERROR: unsupported filter:', filter
            exit(1)

        wig_data[bin_idx] = round(wig_data[bin_idx], PRECISION)


def smooth_dict(d, filter='gaussian', bandwidth=10):
    for key_idx, key in enumerate(d):
        # echo('Smoothing', key)
        smooth_in_place(d[key], filter, bandwidth)

        if key_idx % 10 == 0:
            gc.collect()
        # d[key] = smooth(d[key], filter, bandwidth)


def state_key(s):
    return int(s[1:]) if s[0] == 'E' else s


def read_segmentations(filenames):
    data = {}
    states = set()
    for fname in filenames:
        data[fname] = {}
        echo('Reading:', fname)
        with open_file(fname) as in_f:
            for line in in_f:
                chrom, start, end, state = line.strip().split('\t')
                if chrom not in data[fname]:
                    data[fname][chrom] = []
                data[fname][chrom].append((int(start), int(end), state))
                states.add(state)

    return data, sorted(states, key=state_key)


def filter_chroms(segmentations, chromosomes):
    res = dict((fname, dict()) for fname in segmentations)
    for fname in segmentations:

        for chrom in segmentations[fname]:
            if chrom in chromosomes:
                res[fname][chrom] = segmentations[fname][chrom]
            else:
                echo('WARNING:', fname, 'has a chromosome that other files do not have. skipping:', chrom)
    return res


def normal_pdf(x, mu=0.0, stddev=1.0):
    u = float((x - mu) / stddev)
    y = math.exp(-u * u / 2) / (math.sqrt(2 * math.pi) * stddev)
    return y


def normal_cdf(x, mu=0.0, stddev=1.0):
    return (1 + math.erf((x - mu) / (stddev * math.sqrt(2)))) / 2.


def poisson_pmf(k, Lambda):
    return math.exp(k * math.log(Lambda) - math.lgamma(k + 1.0) - Lambda)


def log_poisson_pmf(k, Lambda):
    return (k * math.log(Lambda) - math.lgamma(k + 1.0) - Lambda) / math.log(2)


def binom_test(observed, n, expected_freq=None, expected=None):

    if expected_freq is None:
        expected_freq = float(expected) / n

    mu = expected_freq * n
    # print mu
    if mu <= 20:
        # use Poisson approximation
        return 1 - sum(poisson_pmf(k, mu) for k in xrange(observed))
    else:
        # use normal approximation
        sigma = math.sqrt(mu * (1 - expected_freq))
        # print sigma
        return 1 - normal_cdf(observed, mu, sigma)


def overlap(s1, e1, s2, e2):
    return min(e1, e2) - max(s1, s2)


def generate_random_chunks(group_A_segmentations, chromosomes, BIN_SIZE, N_RANDOM_CHUNKS=10, _chrom_lengths=None):


    if _chrom_lengths is not None:
        chromosomes = sorted(_chrom_lengths)
        chrom_lengths = [_chrom_lengths[c] for c in chromosomes]

    else:
        chrom_lengths = [group_A_segmentations.values()[0][chrom][-1][1] / BIN_SIZE for chrom in chromosomes]

    total_genome_length = sum(chrom_lengths)

    # echo('Total genome length:', total_genome_length)
    # print chrom_lengths

    RANDOM_CHUNK_LENGTH = 5000000 / BIN_SIZE

    if (N_RANDOM_CHUNKS * RANDOM_CHUNK_LENGTH > total_genome_length / 5 or
        N_RANDOM_CHUNKS * RANDOM_CHUNK_LENGTH > sum(l for l in chrom_lengths if l > RANDOM_CHUNK_LENGTH) / 5):

        return dict((chromosomes[chrom_idx],
                     [(0, chrom_length)]) for chrom_idx, chrom_length in enumerate(chrom_lengths))

    MAX_TRIES = 10000

    restart = True
    random_chunks = None

    while restart:

        restart = False
        random_chunks = {}

        for _ in xrange(N_RANDOM_CHUNKS):
            chunk_start = None
            chrom = None
            tries = 0

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

                    tries += 1
                    if tries > MAX_TRIES:
                        restart = True
                        break

            if restart:
                break

            if chrom not in random_chunks:
                random_chunks[chrom] = []

            random_chunks[chrom].append((chunk_start, chunk_start + RANDOM_CHUNK_LENGTH))
    # print dict((chrom, sorted(random_chunks[chrom])) for chrom in random_chunks)
    return dict((chrom, sorted(random_chunks[chrom])) for chrom in random_chunks)


def slice_segmentations(segmentations_dict, chunk_start, chunk_end):
    return dict((d, segmentations_dict[d][chunk_start:chunk_end]) for d in segmentations_dict)


def chunks(array, chunk_size):
    chunk = []
    for i, el in enumerate(array):
        if i > 0 and i % chunk_size == 0:
            yield chunk
            chunk = []
        chunk.append(el)
    if len(chunk) > 0:
        yield chunk
