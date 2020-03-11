from cyvcf2 import VCF
import sys
from collections import defaultdict
import numpy as np

smp2gt = defaultdict()

vcf = VCF(sys.argv[1], gts012=True)

smp2gt = np.empty(len(vcf.samples), dtype='U10, i4, i4, i4, i4')

for i,smp in enumerate(vcf.samples):
    smp2gt[i][0] = smp

for idx,v in enumerate(vcf('chr17')):
    if idx > 1000: break
    if v.FILTER not in (None, 'PASS'): continue
    if len(v.ALT) > 1: continue
    if len(v.ALT[0]) > 1: continue
    if int(v.INFO.get('AC')) < 2: continue

    gts = v.gt_types

    for i,g in enumerate(gts):
        if g < 0:
            smp2gt[i][-1] += 1
        else:
            smp2gt[i][g+1] += 1

for smp in smp2gt:
    o = 
    print ('\t'.join(map(str, list(smp2gt[smp]))))









