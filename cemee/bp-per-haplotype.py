import sys
import csv
from collections import defaultdict
import glob

chrom = sys.argv[1]

bp_per_hap = defaultdict(int)

for i,f in enumerate(glob.glob('haplotypes-in-lines/indiv-lines/*.chr{}.haplotypes.sorted.csv'.format(chrom))):
    #if i > 0: break
    with open(f, 'r') as fh:
        csvf = csv.reader(fh)
        cur_hap = None
        prev_pos = 0
        for i,l in enumerate(csvf):
            if i == 0: 
                cur_hap = l[3]
                continue
            if l[3] != cur_hap:
                #print (prev_pos, cur_hap, l, bp_per_hap)
                bp_per_hap[cur_hap] += (int(l[0].split(':')[0]) - prev_pos)
                cur_hap = l[3]
                prev_pos = int(l[0].split(':')[0])
print (bp_per_hap)




