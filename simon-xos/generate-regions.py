import argparse
import sys

p = argparse.ArgumentParser()
p.add_argument('--fams')
p.add_argument('-n', type=str, default='1-10')
args = p.parse_args()

with open(args.fams) as fh:
    for i,l in enumerate(fh):
        interval = list(map(int, args.n.split('-')))
        if i + 1 not in range(interval[0], interval[1]): continue
        print >> sys.stderr, "submitting jobs for sample {}".format(l)
        chroms = list(map(str, range(1,23)))
        chroms.append('X')
        for c in chroms:
            print (' '.join([l.rstrip(), c]))
