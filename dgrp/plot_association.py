import sys
from cyvcf2 import VCF
from collections import defaultdict, OrderedDict
import csv
import numpy as np

vcf = VCF(sys.argv[1], gts012=True)

smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))
idx2samp = dict(zip(range(len(vcf.samples)), vcf.samples))

snp_list = []

with open(sys.argv[2], 'r') as fh:
    f = csv.reader(fh, delimiter=' ')
    for i,l in enumerate(f):
        if i == 0: continue
        chrom = l[0].split('_')[0]
        pos = l[0].split('_')[1]
        r, a  = l[1], l[2]
        p = float(l[-1])
        if p > 0.000001: continue
        snp_list.append((chrom, pos, r, a))

phen_dict = defaultdict()
with open(sys.argv[3], 'r') as fh:
    f = csv.reader(fh)
    for l in f:
        phen_dict[l[0]] = l[1]

output = open('snp_distrib.csv', 'a')
for snp in snp_list:

    distribs = defaultdict(list)

    #for v in vcf(snp[0] + ':' + str(int(snp[1]) - 1) + '-' +  snp[1]):
    for v in vcf(snp[0] + ':' + snp[1] + '-' +  snp[1]):
        if len(v.ALT) > 1: continue
        assert v.ALT[0] == snp[2] and v.REF == snp[3]

        for gt in (0,1,2):

            alts = np.where(v.gt_types == gt)[0]

            for idx in alts:
                smp = idx2samp[idx]
                phen = phen_dict[smp.split('_')[-1]]
                distribs[gt].append(phen)

    for gt in distribs:
        counts = distribs[gt]
        for c in counts:
            print (','.join(['_'.join(snp[:2]), str(gt), str(c)]))

    





