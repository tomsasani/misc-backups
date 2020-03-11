from collections import defaultdict
from generate_bst_haps import make_tree
import os

samp2tree = defaultdict()
for f in os.listdir('.'):
    if 'CC' not in f.split('-')[0]: continue
    strain = f.split('-')[0]
    hap_fh = f
    samp2tree[strain] = make_tree(hap_fh)

haps = defaultdict(list)

fh = open('hap_fracs.csv', 'a')
print (','.join(['samp', 'haplotype', 'bp', 'bp_frac']), file=fh)
for s in samp2tree:
    g_size = 0
    per_samp = defaultdict(list)
    tree = samp2tree[s]
    for c in list(map(str, list(range(1,20)))):
        c_t = tree[c].search(0, 250000000)
        for ival in c_t:
            g_size += (ival.end - ival.start)
            h = ival.data
            per_samp[h].append(ival.end - ival.start)

    for h in per_samp:

        total_hap = sum(per_samp[h])
        print (','.join([s, h, str(total_hap), str(float(total_hap) / g_size)]), file=fh)


