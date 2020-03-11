from cyvcf2 import VCF
import argparse
import numpy as np
from pyfaidx import Fasta
from collections import defaultdict
from quicksect import IntervalTree
import itertools
import sys
import time

p = argparse.ArgumentParser()
p.add_argument('--vcf')
p.add_argument('--ref')
p.add_argument('-regions')
p.add_argument('-outname', default='out')
p.add_argument('-exclude')
p.add_argument('-nmer', type=int, default=3)
args = p.parse_args()

def get_regions(path, chrom=None):
    import csv
    regions = defaultdict(list)
    with open(path) as fh:
        f = csv.reader(fh, delimiter='\t')
        for l in f:
            if chrom and l[0] != chrom: continue
            regions[l[0]].append((int(l[1]), int(l[2])))

    return regions

def duplicate_list(l, n):
    return [ele for ele in l for _ in range(n)]

def count_mutations(args, rc, rc_nuc):
    basic_muts = ['C>A','C>G','C>T','A>G','A>C','A>T']

    expanded_muts = []

    if args.nmer == 1:
        expanded_muts = basic_muts
    else:
        k = int((args.nmer - 1) / 2)

        nucleotides = duplicate_list([['A', 'T', 'C', 'G']], k)

        nucleotide_combos = [p for p in itertools.product(*nucleotides)]
        
        for mut in basic_muts:
            for c1 in nucleotide_combos:
                for c2 in nucleotide_combos:
                    c1 = ''.join(c1)
                    c2 = ''.join(c2)
                    expanded_muts.append(c1 + '_' + mut + '_' + c2)

    return expanded_muts

if args.nmer % 2 == 0:
    print ("not a valid nmer count", file=sys.stderr)
    sys.exit()

# read in VCF and REF
vcf = VCF(args.vcf)
reference = Fasta(args.ref)

# read in exclude file and conserved regions
exclude = None
if args.exclude:
    exclude = read_exclude(args.exclude)
region_list = reference.keys()
if args.regions:
    region_list = get_regions(args.regions)

# count numbers of possible mutations
n_mutation_combos = (4 ** (int(args.nmer - 1))) * 6

# get dictionaries of samples and VCF indices
smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))
#smp2idx = {k:v for k,v in smp2idx.items() if k[:2] in ('HG', 'NA')} 
idx2samp = {v:k for k,v in smp2idx.items()}# if k[:2] in ('HG', 'NA')} 

# limit VCF to 1KG samples
vcf = VCF(args.vcf, gts012=True, samples=[s for s in smp2idx])

rc = {'G>T':'C>A', 'G>C':'C>G', 'G>A':'C>T', 'T>C':'A>G', 'T>G':'A>C', 'T>A':'A>T'}
rc_nuc = {'A':'T', 'C':'G', 'T':'A', 'G':'C'}

expanded_muts = count_mutations(args, rc, rc_nuc)

mut2idx = dict(zip(expanded_muts, range(len(expanded_muts))))
idx2mut = {v:k for k,v in mut2idx.items()}
mut_count = 0

ncdng_mutation_counts = np.zeros((len(smp2idx), n_mutation_combos))

r_count = 0

for chrom in region_list:
    v_chrom = chrom
    r_chrom = chrom
    if not any(['chr' in c for c in vcf.seqnames]): v_chrom = chrom.split('chr')[-1]

    r_vcf = vcf(v_chrom)
    for i,v in enumerate(r_vcf):
        if exclude is not None and len(exclude[v.CHROM].search(v.start, v.end)) > 0: continue

        if i % 1000 == 0: print ('done with {} variants'.format(i))
       
        if v.FILTER not in (None, "PASS"): continue
        if len(v.REF) > 1: continue
        if len(v.ALT) > 1: continue
        if len(v.ALT[0]) > 1: continue
        if v.call_rate < 1: continue

        n_nucs = None
        if args.nmer == 1:
            n_nucs = 0
        else: n_nucs = int((args.nmer - 1) / 2)

        if args.nmer == 1:
            upstream, dnstream = None, None
        else:
            upstream = str(reference[r_chrom][v.start - n_nucs:v.start]).upper()
            dnstream = str(reference[r_chrom][v.end:v.end + n_nucs]).upper()

        assert str(reference[r_chrom][v.start]).upper() == v.REF

        if not (len(str(upstream)) == n_nucs and len(str(dnstream)) == n_nucs): continue

        ref, alt = v.REF, v.ALT[0]

        mutation = ref + '>' + alt

        if upstream is not None and ('N' in upstream or 'N' in dnstream): continue

        if mutation in rc and upstream is not None:
            upstream = ''.join([rc_nuc[u] for u in list(upstream)])
            dnstream = ''.join([rc_nuc[d] for d in list(dnstream)])
            mutation = rc[mutation]
            kmer = '_'.join([dnstream, mutation, upstream])
        elif mutation not in rc and upstream is not None:
            kmer = '_'.join([upstream, mutation, dnstream])
        elif mutation in rc and upstream is None:
            mutation = rc[mutation]
            kmer = mutation
        elif mutation not in rc and upstream is None:
            kmer = mutation
        
        kmer_idx = mut2idx[kmer]

        gts = v.gt_types

        known_gts = np.isin(gts, [0,1,2]) 
        known_gt_idx = np.where(known_gts)[0]
        ncdng_mutation_counts[known_gt_idx[:,None],kmer_idx] += gts[known_gt_idx].reshape(known_gt_idx.shape[0], 1)
        
mlist = ','.join(mut2idx.keys())

fh = 'k{}.{}.mutation_spectra.txt'.format(args.nmer, args.outname)

fh = open(fh, mode='a')
print (','.join(['samp', mlist]), file=fh)
for i,samp in enumerate(ncdng_mutation_counts):
    n_vals = list(map(int, list(samp)))
    n_vals = ','.join(list(map(str, n_vals)))
    print (','.join([idx2samp[i], n_vals]), file=fh)
