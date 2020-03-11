from cyvcf2 import VCF
import argparse
import numpy as np
from pyfaidx import Fasta
from collections import defaultdict
import itertools
import sys

p = argparse.ArgumentParser()
p.add_argument('--vcf')
p.add_argument('--ref')
p.add_argument('-samples')
p.add_argument('-nmer', type=int, default=3)
p.add_argument('-r')
args = p.parse_args()

if args.nmer % 2 == 0:
    print ("not a valid nmer count", file=sys.stderr)
    sys.exit()

if args.samples:
    vcf = VCF(args.vcf, samples=args.samples.split(','))
else:
    vcf = VCF(args.vcf)

smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))

print ('{} samples in the VCF'.format(len(smp2idx)), file=sys.stderr)

HOM_ALT = 3
HET = 1
HOM_REF = 0
UNK = 2

reference = Fasta(args.ref)

n_mutation_combos = (4 ** (int(args.nmer - 1))) * 6

mutation_counts = {s:np.zeros(n_mutation_combos) for s in vcf.samples}

basic_muts = ['C>A','C>G','C>T','A>G','A>C','A>T']

rc = {'G>T':'C>A', 'G>C':'C>G', 'G>A':'C>T', 'T>C':'A>G', 'T>G':'A>C', 'T>A':'A>T'}

rc_nuc = {'A':'T', 'C':'G', 'T':'A', 'G':'C'}

def duplicate_list(l, n):
    return [ele for ele in l for _ in range(n)]

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

mut2idx = dict(zip(expanded_muts, range(len(expanded_muts))))
mut_count = 0

if args.r:
    vcf = vcf(args.r)

for i,v in enumerate(vcf):
    if i % 10000 == 0: print ('done with {} variants, on chrom {}'.format(i, v.CHROM), file=sys.stderr)

    #if v.FILTER not in (None, "PASS"): continue

    if len(v.REF) > 1: continue
    if len(v.ALT) > 1: continue
    if len(v.ALT[0]) > 1: continue
    if v.INFO.get('MQ') is not None and int(v.INFO.get('MQ')) < 30: continue
    if v.CHROM in ("chrX", "chrY", "X", "Y"): continue
    if v.CHROM not in ['chr' + str(n) for n in range(1,20)]: continue

    #n_no_depth = np.where(v.format('AD') < 0)[0].shape[0]

    ### TODO: why does this filter out so many? why do so many have GT but no DP? imputed?
    #if float(n_no_depth) / (len(smp2idx) * 2) > 0.05: continue    

    gts = v.gt_types

    call_rate = np.sum(gts != 2) / float(len(smp2idx))

    if call_rate < 0.99: continue

    n_nucs = None
    if args.nmer == 1:
        n_nucs = 0
    else: n_nucs = int((args.nmer - 1) / 2)

    if args.nmer == 1:
        upstream, dnstream = None, None
    else:
        upstream = str(reference[v.CHROM][v.start - n_nucs:v.start]).upper()
        dnstream = str(reference[v.CHROM][v.start:v.start + n_nucs]).upper()

    assert str(reference[v.CHROM][v.start]).upper() == v.REF

    ref, alt = v.REF, v.ALT[0]

    mutation = ref + '>' + alt

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
    
    for samp in smp2idx:
        gt = gts[smp2idx[samp]]
        if gt == HOM_ALT: 
            mutation_counts[samp][kmer_idx] += 2
        elif gt == HET:
            mutation_counts[samp][kmer_idx] += 1
        else: 
            mut_count += 1
    #print (np.where(mutation_counts['ECA744']>0)[0].shape)

mlist = ','.join(mut2idx.keys())

print (','.join(['samp', mlist]))
for samp in mutation_counts:
    vals = list(map(int, list(mutation_counts[samp])))
    vals = ','.join(list(map(str, vals)))
    print (','.join([samp, vals]))

print ('total of {} mutations that pass'.format(mut_count), file=sys.stderr)
