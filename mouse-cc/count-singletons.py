from cyvcf2 import VCF
import sys
import numpy as np
import gzip
import csv
import glob
import os
import fnmatch
from collections import defaultdict, Counter
from generate_bst_haps import make_tree

vcf = VCF(sys.argv[1], gts012=True)

samps = list(vcf.samples)

s_counts_snp = np.zeros((len(vcf.samples), 36))
s_counts_indel = np.zeros((len(vcf.samples), 36))

idx2samp = dict(zip(range(len(vcf.samples)), vcf.samples))
samp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))

hap2idx = defaultdict(int)
haps_counted = 0

count_per_lines = defaultdict(int)
callable_per_lines = np.zeros(len(vcf.samples), dtype=int)
unk_per_lines = np.zeros(len(vcf.samples), dtype=int)

samp2tree = defaultdict()

for s in list(vcf.samples):
    hap_fh = None
    for f in os.listdir('.'):
        if fnmatch.fnmatch(f, s + '*'):
            hap_fh = f
    samp2tree[s] = make_tree(hap_fh)

ii = np.arange(len(vcf.samples))
### TODO: make sure there are intervals for the ends of chromosomes!
for i,v in enumerate(vcf):
    if i % 10000 == 0: print ('done with {} variants, on {}:{}'.format(i, v.CHROM, v.POS))
#    if v.var_type != "snp": continue

    gts = v.gt_types
    unk_gts = np.where(gts == 3)

    ad, rd = v.gt_alt_depths, v.gt_ref_depths
    gq = v.gt_quals
    td = ad + rd

    good_gq = np.where(gq > 20)
    good_depth = np.where(td >= 15)

    good_sites = np.intersect1d(good_gq, good_depth)
    
    if good_sites.shape[0] == 0: continue

    callable_per_lines[np.intersect1d(good_sites, np.where(gts != 3))] += 1

    if v.INFO.get('AC') != 2: continue
    
    unique, counts = np.unique(gts, return_counts=True)
    uc = dict(zip(unique, counts))
    if 2 not in uc: continue
    
    singleton = np.where(gts == 2)[0][0]
    
    if not gq[singleton] > 20: continue
    if not td[singleton] >= 15: continue
    if not (ad[singleton] / float(td[singleton])) == 1.: continue

    cc_line = idx2samp[singleton]

    hap_tree = samp2tree[cc_line]

    founder_hap = hap_tree[str(v.CHROM).split('chr')[-1]].search(v.start, v.end)

    if len(founder_hap) == 0: continue
    founder_hap = str(founder_hap[0].data)

    # make sure that there's at least one other mouse
    # with the same founder haplotype, but not the ALT

    # NOTE: currently skipping this, since we're assuming
    # all variants in `Private_Variants.csv` are "good"
    
    samp_with_hap = 0
    samp_with_hap_hr = 0
    for mouse in samp2tree:
        if mouse == cc_line: continue
        f_hap = samp2tree[mouse][str(v.CHROM).split('chr')[-1]].search(v.start, v.end)
        if len(f_hap) == 0: continue
        f_hap = str(f_hap[0].data)

        if f_hap != founder_hap: continue
        samp_with_hap += 1
        m_idx = samp2idx[mouse]
        if not gts[m_idx] == 0: continue
        if not (gq[m_idx] >= 20 and td[m_idx] >= 15): continue
        samp_with_hap_hr += 1

    print ('{} mice with same haplotype, {} that are HQ HOM_REF'.format(samp_with_hap, samp_with_hap_hr))
    
    # make sure that at least one other sample has the same 
    # haplotype at this site

    if samp_with_hap_hr < 1: continue
    

    if founder_hap is None: continue
    founder_hap_idx = None
    if founder_hap in hap2idx:
        founder_hap_idx = hap2idx[founder_hap]
    else:
        hap2idx[founder_hap] = haps_counted
        founder_hap_idx = haps_counted
        haps_counted += 1

    if v.var_type == 'snp':
        s_counts_snp[singleton][founder_hap_idx] += 1
    else:
        s_counts_indel[singleton][founder_hap_idx] += 1


kk = [*hap2idx]
kk_add = ['unk' + str(i) for i in range(36 - len(kk))]
kk.extend(kk_add)
kk.insert(0, 'samp')

def write_singletons(array, fname, idx2samp, header=[]):
    fh = open(fname, mode='a')
    print (','.join(header), file=fh)
    print (np.sum(array))
    for i,ril in enumerate(array):
        counts = list(map(int, list(ril)))
        counts = ','.join(list(map(str, counts)))
        print (','.join([idx2samp[i], counts]), file=fh)

write_singletons(s_counts_snp, 'singleton_counts_snp.csv', idx2samp, header=kk)
write_singletons(s_counts_indel, 'singleton_counts_indel.csv', idx2samp, header=kk)
