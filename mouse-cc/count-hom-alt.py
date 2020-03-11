import sys
from cyvcf2 import VCF

vcf = VCF(sys.argv[1])

smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))


count = 0 
for v in vcf:
    if v.var_type != "snp": continue
    if v.call_rate < 1: continue
    gts = v.gt_types

    if gts[smp2idx['CC010']] != 3: continue
    if 'Low' in v.FILTER: continue

    print (str(v).rstrip())
    count += 1

print (count)


