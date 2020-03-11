import sys
from cyvcf2 import VCF
import numpy as np
from collections import defaultdict

vcf = VCF(sys.argv[1])

idx2samp = dict(zip(range(len(vcf.samples)), vcf.samples))

unk_counts = np.zeros(len(idx2samp), dtype=np.float64)
all_counts = np.zeros(len(idx2samp), dtype=np.float64)

for i,v in enumerate(vcf('chr17:1-1000000')):
    if i % 10000 == 0: print ('done with {} variants'.format(i), file=sys.stderr)

    if v.FILTER not in (None, 'PASS'): continue
    if len(v.ALT) > 1: continue
    if len(v.ALT[0]) > 1: continue
    if len(v.REF) > 1: continue

    gts = v.gt_types

    unk_idx = np.where(gts == 2)[0]

    unk_counts[unk_idx] += 1
    all_counts[np.arange(0,len(idx2samp))] += 1

unk_frac = unk_counts / all_counts

all_unk = np.where(unk_frac == 1)[0]

unk_samps = '\n'.join([idx2samp[idx] for idx in idx2samp if idx in list(all_unk)])

print (unk_samps)
