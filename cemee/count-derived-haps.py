from cyvcf2 import VCF
import sys
import numpy as np
import gzip
import csv
import glob
import os
from collections import defaultdict, Counter

vcf = VCF(sys.argv[1])
s_counts = np.zeros((len(vcf.samples), 16))
idx2samp = dict(zip(range(len(vcf.samples)), vcf.samples))
samp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))

hap2idx = defaultdict(int)
haps_counted = 0

def bp_per_hap(chrom):
    bp_per_hap = defaultdict(lambda: defaultdict(int))
    for i,f in enumerate(glob.glob('haplotypes-in-lines/indiv-lines/*.chr{}.haplotypes.sorted.csv'.format(chrom))):
        line = f.split('/')[-1].split('.')[0]
        with open(f, 'r') as fh:
            csvf = csv.reader(fh)
            cur_hap = None
            prev_pos = 0
            for i,l in enumerate(csvf):
                if i == 0: 
                    cur_hap = l[3]
                    continue
                if l[3] != cur_hap:
                    bp_per_hap[line][cur_hap] += (int(l[0].split(':')[0]) - prev_pos)
                    cur_hap = l[3]
                    prev_pos = int(l[0].split(':')[0])
    return bp_per_hap

def get_hap(founder, pos, chrom='1'):
    f = 'haplotypes-in-lines/indiv-lines/' + founder + '.' + '{}'.format(chrom) + '.haplotypes.sorted.csv'
    if not os.path.exists(f): return None
    with open(f, 'rt') as ff:
        fh = csv.reader(ff)
        cur_hap = None
        for i,l in enumerate(fh):
            cur_pos = int(l[0].split(':')[0])
            cur_hap = l[3]
            if cur_pos < pos: 
                cur_hap = l[3]
                continue
            else: return cur_hap

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

s2l = sra2line(sys.argv[2])

chromd = {'I':'chr1', 'II':'chr2', 'III':'chr3', 'IV':'chr4', 'V':'chr5', 'X':'chrX'}

count_per_lines = defaultdict(int)
callable_per_lines = np.zeros(len(vcf.samples))

ii = np.arange(len(vcf.samples))

hap_counts = defaultdict(int)
for i,v in enumerate(vcf('I')):
    if i % 10000 == 0: print ('done with {} variants, on chr{}:{}'.format(i, v.CHROM, v.POS))

    if v.FILTER not in ('PASS', None): continue
    if v.var_type != "snp": continue
    
    gts = v.gt_types
    ad, rd = v.gt_alt_depths, v.gt_ref_depths
    gq = v.gt_quals
    td = ad + rd

    good_gq = np.where(gq >= 20)
    good_depth = np.where(td >= 10)

    good_sites = np.intersect1d(good_gq, good_depth)
    good_idx = np.intersect1d(ii, good_sites)

    if good_idx.shape == 0: continue
    callable_per_lines[good_idx] += 1

    if v.INFO.get('AC') > 2: continue
    if v.INFO.get('AC') == 0: continue
    unique, counts = np.unique(gts, return_counts=True)
    uc = dict(zip(unique, counts))
    if 3 in uc and 1 in uc: continue
    if 1 in uc and uc[1] > 1: continue
    
    singleton = np.where((gts == 1) | (gts == 3))[0][0]
    if not gq[singleton] >= 20: continue
    if not td[singleton] >= 10: continue
    sra_line = idx2samp[singleton]
    
    founder = s2l[sra_line]
    try:
        line_num = founder.split('L')[1]
    except IndexError: continue

    founder_hap = get_hap(founder, int(v.POS), chrom=chromd[v.CHROM])
    
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
    
    count_per_lines[sra_line] += 1

bp_per_hap = bp_per_hap('1')

idx2hap = {v:k for k,v in hap2idx.items()}

fh = open('haps_vs_singletons.csv', mode='a')
kk = [*hap2idx]
kk.insert(0, 'line')
print (','.join(kk), file=fh)
for i,ril in enumerate(s_counts):
    counts = list(map(int, list(ril)))
    adj_counts = []
    for ii,c in enumerate(counts):
        if ii not in idx2hap:
            adj_counts.append(0)
            continue
        if c == 0:
            adj_counts.append(0)
            continue
        orig_founder_hap = idx2hap[ii]
        sra_name = idx2samp[i]
        ril_name = s2l[sra_name]
        orig_founder_bp = bp_per_hap[ril_name][orig_founder_hap]
        print (orig_founder_hap, orig_founder_bp, c)
        if c > 0 and orig_founder_bp == 0:
            adj_counts.append(0)
            continue
        adj_counts.append(c / float(orig_founder_bp))
    counts = ','.join(list(map(str, counts)))
    adj_counts = ','.join(list(map(str, adj_counts)))
    print (','.join([s2l[idx2samp[i]] + '_adj', adj_counts]), file=fh)
    print (','.join([s2l[idx2samp[i]], counts]), file=fh)
