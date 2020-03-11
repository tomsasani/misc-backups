from cyvcf2 import VCF
import sys
import numpy as np
import gzip
import csv
import glob
import os
from collections import defaultdict, Counter

def get_hap(founder, pos, chrom='1'):
    f = 'haplotypes-in-lines/indiv-lines/' + founder + '.' + 'chr{}'.format(chrom) + '.haplotypes.sorted.csv'
    if not os.path.exists(f): return None
    with open(f, 'rt') as ff:
        fh = csv.reader(ff)
        hap = None
        for i,l in enumerate(fh):
            cur_pos = int(l[0].split(':')[0])
            if cur_pos != pos: continue
            hap = l[3]
            return hap

def sra2line(filename):
    s2l = defaultdict(str)
    with open(filename, 'r') as fh:
        f = csv.reader(fh)
        for i,l in enumerate(f):
            if i == 0: continue
            if len(l) == 0: continue
            sra = l[0]
            line = l[29]
            s2l[sra] = line
    return s2l

vcf = VCF(sys.argv[1], gts012=True)
s_counts = np.zeros((len(vcf.samples), 16))
idx2samp = dict(zip(range(len(vcf.samples)), vcf.samples))
samp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))

hap2idx = defaultdict(int)
haps_counted = 0

#s2l = sra2line(sys.argv[2])

chromd = {'I':'chr1', 'II':'chr2', 'III':'chr3', 'IV':'chr4', 'V':'chr5', 'X':'chrX'}

count_per_lines = defaultdict(int)
callable_per_lines = np.zeros(len(vcf.samples), dtype=int)
unk_per_lines = np.zeros(len(vcf.samples), dtype=int)

ii = np.arange(len(vcf.samples))

hap_counts = defaultdict(int)
for i,v in enumerate(vcf):
    if i % 10000 == 0: print ('done with {} variants, on chr{}:{}'.format(i, v.CHROM, v.POS))
    #if v.FILTER not in ('PASS', None): continue
    if v.call_rate < 1: continue
    if v.var_type != "snp": continue

    gts = v.gt_types
    unk_gts = np.where(gts == 3)
    unk_per_lines[unk_gts] += 1

    if v.INFO.get('AC') != '2': continue
    unique, counts = np.unique(gts, return_counts=True)
    uc = dict(zip(unique, counts))
    if 2 not in uc: continue
    
    singleton = np.where(gts == 2)[0][0]
    
    founder = idx2samp[singleton]
    
    try:
        line_num = founder.split('L')[1]
    except IndexError: continue
    founder_hap = get_hap(founder, int(v.POS), chrom=v.CHROM)
    if founder_hap is None: continue
    founder_hap_idx = None
    if founder_hap in hap2idx:
        founder_hap_idx = hap2idx[founder_hap]
    else:
        hap2idx[founder_hap] = haps_counted
        founder_hap_idx = haps_counted
        haps_counted += 1
    s_counts[singleton][founder_hap_idx] += 1
    hap_counts[founder_hap] += 1
    if s_counts.sum() % 100. == 0:
        print ('counted {} singletons'.format(s_counts.sum()))
        print (Counter(hap_counts))
    
    count_per_lines[founder] += 1

fh = open('haps_vs_singletons.imputed.csv', mode='a')
kk = [*hap2idx]
kk.insert(0, 'sra')
print (','.join(kk), file=fh)
for i,ril in enumerate(s_counts):
    counts = list(map(int, list(ril)))
    counts = ','.join(list(map(str, counts)))
    print (','.join([idx2samp[i], counts]), file=fh)


### need to do as function of callable sites in each strain
fh = open('singletons_in_lines.imputed.csv', mode='a')
print (','.join(['sra', 'line', 'count', 'callable_sites', 'unk_sites']), file=fh)
for sra_l in count_per_lines:
    count = count_per_lines[sra_l]
    line = sra_l
    idx = samp2idx[sra_l]
    print (','.join([sra_l, line, str(count), str(callable_per_lines[idx]), str(unk_per_lines[idx])]), file=fh)



