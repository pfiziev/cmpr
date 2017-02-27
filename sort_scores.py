import argparse
import sys
from utils import open_file, echo

parser = argparse.ArgumentParser()

parser.add_argument('-q', nargs='+', help='q-value scores')

args = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    exit(1)

for fname in args.q:
    out_fname = fname.replace('.gz', '').replace('.bed', '.sorted_by_score.bed.gz')
    echo('Processing:', fname)

    with open_file(fname) as in_f, \
         open_file(out_fname, 'w') as out_f:

        for _, line in sorted([(float(line.split()[4]), line) for line in in_f], reverse=True):
            out_f.write(line)
