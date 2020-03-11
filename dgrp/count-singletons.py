from cyvcf2 import VCF
import sys
import numpy as np
import gzip
import csv
import glob
import os
import fnmatch
from collections import defaultdict, Counter

vcf = VCF(sys.argv[1], gts012=True)

samps = list(vcf.samples)

s_counts_snp = np.zeros((len(vcf.samples)))
s_counts_indel = np.zeros((len(vcf.samples)))

idx2samp = dict(zip(range(len(vcf.samples)), vcf.samples))
samp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))

ii = np.arange(len(vcf.samples))
### TODO: make sure there are intervals for the ends of chromosomes!
for i,v in enumerate(vcf):
    if i % 10000 == 0: 
        print ('done with {} variants, on {}:{}'.format(i, v.CHROM, v.POS))
        print (np.sum(s_counts_snp))
    if v.call_rate != 1.: continue
    gts = v.gt_types
    unk_gts = np.where(gts == 3)

    unique, counts = np.unique(gts, return_counts=True)
    uc = dict(zip(unique, counts))
    if 2 not in uc: continue
    
    singleton = np.where(gts == 2)[0][0]
    
    if v.var_type == 'snp':
        s_counts_snp[singleton] += 1
    else:
        s_counts_indel[singleton] += 1

def write_singletons(array, fname, idx2samp, header=[]):
    fh = open(fname, mode='a')
    print (np.sum(array))
    for i,ct in enumerate(array):
        print (','.join([idx2samp[i].split('_')[-1], str(int(ct))]), file=fh)

write_singletons(s_counts_snp, 'singleton_counts_snp.csv', idx2samp)
write_singletons(s_counts_indel, 'singleton_counts_indel.csv', idx2samp)
