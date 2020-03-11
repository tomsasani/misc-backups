from quicksect import IntervalTree
from collections import defaultdict
import sys
import gzip
import toolshed as ts
import numpy as np

def make_tree(path):
    tree = defaultdict(IntervalTree)
    prev_chrom, prev_pos = None, 0
    prev_hap = None
    added = defaultdict(int)
    for i, line in enumerate(ts.reader(path, sep=',')):
        chrom = line['chromosome']
        pos = line['position(B38)']
        hap_probs = list(line.items())[3:]
        hap_probs_np = np.array([v for k,v in hap_probs])
        max_hap = [v for i,v in enumerate(hap_probs) if i == np.argmax(hap_probs_np)][0]
        hap, score = max_hap
        if float(score) < 0.8: continue
        if i == 0:
            prev_hap = hap
            prev_chrom = chrom
            prev_pos = pos
        if chrom == prev_chrom and hap != prev_hap:
            tree[prev_chrom].add(int(prev_pos), int(pos), other=prev_hap)
#            print ('added {}:{}-{}'.format(prev_chrom, prev_pos, pos))
            prev_pos = pos
            prev_chrom = chrom
            prev_hap = hap
            added[chrom] = int(pos)
        elif chrom != prev_chrom:
            tree[prev_chrom].add(int(prev_pos), added[prev_chrom], other=prev_hap)
#            print ('added {}:{}-{}'.format(prev_chrom, prev_pos, added[prev_chrom]))
            prev_pos = 0
            prev_chrom = chrom
            prev_hap = hap
        else: 
            added[chrom] = int(pos)
            continue
    return tree


def run():
    import pickle
    tree = make_tree(sys.argv[1])

    print (tree['X'].search(0,180000000))

if __name__ == "__main__":
    run()


