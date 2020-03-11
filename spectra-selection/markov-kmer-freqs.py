import numpy as np
from pyfaidx import Fasta
import argparse
from collections import defaultdict, Counter
from quicksect import IntervalTree
import scipy.stats as ss
import csv
import sys
import re
import gzip

p = argparse.ArgumentParser()
p.add_argument('--ref')
p.add_argument('--sequence')
p.add_argument('-regions')
p.add_argument('-pctile')
p.add_argument('-chrom', default='19')
p.add_argument('-n', default=6, type=int)
args = p.parse_args()

reference = Fasta(args.ref)

chroms = [] 
for c in reference.keys():
    chroms.append(c)

chrom = args.chrom
if len(chroms) == 1:
    chrom = chroms[0]

def read_exclude(path, chrom=None):
    if path is None:
        return None
    tree = defaultdict(IntervalTree)
    # in case they sent a region.
    if chrom is not None:
        chrom = chrom.split(":")[0]

    added = 0
    for i, line in enumerate((gzip.open if path.endswith(".gz") else open)(path)):
        toks = line.rstrip().split("\t")
        # skip header if necessary.
        if i == 0:
            try:
                int(toks[1])
            except ValueError:
                continue
        if chrom is not None and toks[0] != chrom:
            continue
        added += 1
        tree[toks[0]].add(int(toks[1]), int(toks[2]))
    if added == 0:
        sys.stderr.write("didnt add any intervals to exclude for chrom:%s\n" % chrom)
    return tree

def get_coords(my_f, dd, pctile):
    with open(my_f) as f:
        fh = csv.reader(f, delimiter='\t')
        for l in fh:
            if float(l[-1]) < pctile: continue
            dd[l[0]].append([int(l[1]), int(l[2])])
    return dd

def count_tfbs(f, sequence, chrom=None, regions=None):

    w_size = len(sequence)
    increm = 1

    seq = None
    
    seq = str(f[chrom])
    
    rr = None

    if regions is not None:
        rr = regions[chrom]
    else:
        rr = [(0,len(str(seq)))]

    tata_counts = 0
    total_seq = 0

    for r in rr:
        i = 0
        if r[1] - r[0] == 0: continue
        r_seq = seq[r[0]:r[1]]
        if len(r_seq) == 0: continue
        
        while i + w_size <= len(r_seq): 
            word = str(r_seq[i:i+w_size]).upper()
            tata_in_word = re.search(sequence, word)
            if tata_in_word is not None: 
                tata_counts += 1
            i += 1
        total_seq += (r[1] - r[0])

    return tata_counts, total_seq

def exp_mm(f, sequence, chrom=None, regions=None):


    import exrex
    seqs_of_interest = list(exrex.generate(sequence))
    nuc_freq = [len(s) for s in seqs_of_interest][0] - 1

    seq = str(f[chrom])

    rr = None
    if regions is not None:
        rr = regions[chrom]
    else:
        rr = [(0,len(str(seq)))]
    

    print ('considering {}'.format(','.join(seqs_of_interest)))

    all_nmers = defaultdict(int)
    all_sub_nmers = defaultdict(int)

    for r in rr:
        r_seq = seq[r[0]:r[1]]
        if len(r_seq) == 0: continue

        i = 0
        while i + nuc_freq <= len(r_seq):
            word = str(r_seq[i:i+nuc_freq]).upper()
            all_nmers[word] += 1
            i += 1

        i = 0
        while i + (nuc_freq - 1) <= len(r_seq):
            s = r_seq[i:i+nuc_freq - 1]
            all_sub_nmers[s] += 1
            i += 1

    exp_seq_frequencies = defaultdict(float)

    for seqoi in seqs_of_interest:

        nmers = []
        sub_nmers = []

        i = 0
        while i + nuc_freq <= len(seqoi):
            nmers.append(seqoi[i:i+nuc_freq])
            i += 1
      
        i = 1 
        while i > 0 and (i + nuc_freq - 1) < len(seqoi):
            sub_nmers.append(seqoi[i:i+nuc_freq - 1])
            i += 1

        num, denom = None, None

        for nmer in nmers:
            f = all_nmers[nmer]
            if num is None:
                num = f
            else: num = num * f

        for sub_nmer in sub_nmers:
            f = all_sub_nmers[sub_nmer]
            if denom is None:
                denom = f
            else: denom = denom * f

        exp_seq_frequencies[seqoi] = num/denom
    return exp_seq_frequencies

regions = None
if args.regions:
    regions = defaultdict(list)
    regions = get_coords(args.regions, regions, float(args.pctile))

exp_tata = exp_mm(reference, args.sequence, chrom=chrom, regions=regions)
summed_exp = sum([f[1] for f in exp_tata.items()])

obs_tata, total_seq = count_tfbs(reference, args.sequence, chrom=chrom, regions=regions) 

g, p, dof, expctd = ss.chi2_contingency([[obs_tata, summed_exp],[total_seq, total_seq]], lambda_='log-likelihood')

enrich = 'n/a'

if obs_tata / summed_exp < 1: enrich = 'UNDER'
elif obs_tata / summed_exp > 1: enrich = 'OVER'

print (','.join(map(str, [chrom, obs_tata, int(summed_exp), total_seq, enrich, p])))
