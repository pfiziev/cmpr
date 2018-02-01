import gzip
import sys
import io
import datetime

import math


def read_p_values(fnames, n_states=None):
    data = {}
    for fname in fnames:
        print >>sys.stderr, datetime.datetime.now(), 'Reading:', fname

        # with io.BufferedReader(gzip.open(fname)) as in_f:
        with open(fname) as in_f:
            for l in in_f:
                buf = l.strip().split()
                key = '\t'.join(buf[:3])
                value = float(buf[8])
                # data[key] = value
                if key not in data:
                    data[key] = 1

                if n_states is not None:
                    data[key] = min(min(1, value * n_states), data[key])
                else:
                    data[key] = min(value, data[key])

    return data


if len(sys.argv) == 1:
    print 'usage: %s n_states output.png overall-score-file per-state-score-files' % __file__
    exit(1)

n_states = int(sys.argv[1])
out_fname = sys.argv[2]
overall_fname = sys.argv[3]
per_state_fnames = sys.argv[4:]

overall_p_values = read_p_values([overall_fname])
per_state_p_values = read_p_values(per_state_fnames, n_states=n_states)

if len(overall_p_values) != len(per_state_p_values):
    print 'error - lengths do not match:', len(overall_p_values), len(per_state_p_values)

keys = sorted(overall_p_values)

# for key in keys:
#     if overall_p_values[key] < per_state_p_values[key]:
#         print '\t'.join(map(str, [key, overall_p_values[key], per_state_p_values[key]]))

for key in keys:
    print '\t'.join(map(str, [key, overall_p_values[key], per_state_p_values[key]]))

# exit(1)

print >>sys.stderr, datetime.datetime.now(), 'plotting'

import matplotlib

matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
plt.figure(figsize=(6, 5))
# plt.title(overall_fname)
min_p = min(min(v for v in overall_p_values.values() if v > 0),
            min(v for v in per_state_p_values.values() if v > 0))
min_p /= 2

log10 = lambda v: -math.log(v if v > 0 else min_p, 10)
# log10 = lambda v: v #-math.log(v if v > 0 else min_p, 10)

x = [log10(overall_p_values[k]) for k in keys]
y = [log10(per_state_p_values[k]) for k in keys]

plt.hist2d(x, y, bins=100, norm=LogNorm())

l1 = min(min(x), min(y))
l2 = max(max(x), max(y))
plt.xlim((l1, l2))
plt.ylim((l1, l2))


plt.xlabel('overall q-values')
plt.ylabel('per-state q-values')

plt.plot([min(x), max(x)], [min(y), max(y)], color='black')
plt.colorbar()
plt.savefig(out_fname)
