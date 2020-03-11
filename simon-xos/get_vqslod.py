from cyvcf2 import VCF
import sys

vqslod = []

for i,v in enumerate(VCF(sys.argv[1])):
    if v.FILTER not in ('PASS', None): continue
    if v.call_rate < 0.9: continue
    if len(v.ALT) > 1: continue
    if len(v.ALT[0]) > 1: continue
    if len(v.REF) > 1: continue

    vqslod.append(float(v.INFO.get('VQSLOD')))


    if i > 100000: break

import numpy as np
print (np.histogram(vqslod))


