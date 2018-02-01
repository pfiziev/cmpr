import os
import sys
from scipy.stats import pearsonr
import datetime

from utils import open_file
from pybedtools import *

BIN_SIZE = 200

if len(sys.argv) == 1:
    print 'usage: %s target-feature.bed score-files.bed' % __file__
    exit(1)

target_bed = BedTool(sys.argv[1])

for fname in sys.argv[2:]:
    score = []
    delta_expression = []
    for f in target_bed.intersect(BedTool(fname), wo=True):
        score.append(float(f[4]))
        delta_expression.append(float(f[10]))

    print datetime.datetime.now(), os.path.split(sys.argv[1])[1], fname, pearsonr(score, delta_expression)[0]

