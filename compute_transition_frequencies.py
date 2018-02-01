import sys
from utils import *

if len(sys.argv) == 1:
    print 'usage: %s silvermark-summary.bed.gz q-value-threshold' % __file__
    exit(1)

qvalue = float(sys.argv[2])
transitions_A_to_B = {}

states = set()

with open_file(sys.argv[1]) as in_f:
    for l in in_f:
        buf = l.strip().split()
        if float(buf[8]) > qvalue:
            break

        s1, s2 = buf[6].split('-')
        if s1 not in transitions_A_to_B:
            transitions_A_to_B[s1] = {}
        if s2 not in transitions_A_to_B[s1]:
            transitions_A_to_B[s1][s2] = 0

        transitions_A_to_B[s1][s2] += 1

        states.add(s1)
        states.add(s2)

states = sorted(states, key=state_key)
print 'A -> B'
print '\t'.join(['state'] + states)
for s1 in sorted(transitions_A_to_B, key=state_key):
    total = float(sum(transitions_A_to_B[s1].values()))
    print '\t'.join([s1] + ['%.2lf' % (100 * transitions_A_to_B[s1].get(s2, 0) / total) for s2 in states])

print
print 'B -> A'
states = sorted(states, key=state_key)
print '\t'.join(['state'] + states)
for s1 in states:
    total = float(sum(transitions_A_to_B[s].get(s1, 0) for s in states))
    print '\t'.join([s1] + ['%.2lf' % (100 * transitions_A_to_B[s2].get(s1, 0) / total) for s2 in states])

print
print 'counts'
states = sorted(states, key=state_key)
print '\t'.join(['state'] + states)
for s1 in states:
    #total = float(sum(transitions_A_to_B[s].get(s1, 0) for s in states))
    print '\t'.join([s1] + ['%d' % transitions_A_to_B[s1].get(s2, 0) for s2 in states])


