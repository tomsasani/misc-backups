from cyvcf2 import VCF
import argparse
import numpy as np
from pyfaidx import Fasta
from collections import defaultdict
import sys

p = argparse.ArgumentParser()
p.add_argument('--vcf')
p.add_argument('--ref')
p.add_argument('-samples')
args = p.parse_args()

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

mutation_counts = {s:np.zeros(96) for s in vcf.samples}

basic_muts = ['C>A','C>G','C>T','A>G','A>C','A>T']

rc = {'G>T':'C>A', 'G>C':'C>G', 'G>A':'C>T', 'T>C':'A>G', 'T>G':'A>C', 'T>A':'A>T'}

rc_nuc = {'A':'T', 'C':'G', 'T':'A', 'G':'C'}

expanded_muts = []

for mut in basic_muts:
    for nuc1 in ['A','T','C','G']:
        for nuc2 in ['A','T','C','G']:
            expanded_muts.append(nuc1 + '_' + mut + '_' + nuc2)

mut2idx = dict(zip(expanded_muts, range(len(expanded_muts))))
mut_count = 0
for i,v in enumerate(vcf):
    if i % 10000 == 0: print ('done with {} variants, on chrom {}'.format(i, v.CHROM), file=sys.stderr)
    if v.FILTER not in (None, "PASS"): continue

    if len(v.REF) > 1: continue
    if len(v.ALT) > 1: continue
    if len(v.ALT[0]) > 1: continue
    if v.INFO.get('MQ') is not None and int(v.INFO.get('MQ')) < 30: continue

    #n_no_depth = np.where(v.format('AD') < 0)[0].shape[0]

    ### TODO: why does this filter out so many? why do so many have GT but no DP? imputed?
    #if float(n_no_depth) / (len(smp2idx) * 2) > 0.05: continue    

    gts = v.gt_types

    call_rate = np.sum(gts != 2) / float(len(smp2idx))

    if call_rate < 0.99: continue

    upstream = str(reference[v.CHROM][v.start - 1]).upper()
    dnstream = str(reference[v.CHROM][v.start + 1]).upper()

    assert str(reference[v.CHROM][v.start]).upper() == v.REF

    ref, alt = v.REF, v.ALT[0]

    mutation = ref + '>' + alt

    three_mer = '_'.join([upstream, mutation, dnstream])

    if mutation in rc:
        upstream = rc_nuc[upstream]
        dnstream = rc_nuc[dnstream]
        mutation = rc[mutation]
        three_mer = '_'.join([dnstream, mutation, upstream])

    three_mer_idx = [v for k,v in mut2idx.items() if k == three_mer][0]

    
    for samp in smp2idx:
        gt = gts[smp2idx[samp]]
        if gt == HOM_ALT: 
            mutation_counts[samp][three_mer_idx] += 2
        elif gt == HET:
            mutation_counts[samp][three_mer_idx] += 1
        else: 
            mut_count += 1

mlist = ','.join(mut2idx.keys())

print (','.join(['samp', mlist]))
for samp in mutation_counts:
    vals = list(map(int, list(mutation_counts[samp])))
#    vals = [v / float(sum(vals)) for v in vals]
    vals = ','.join(list(map(str, vals)))
    print (','.join([samp, vals]))

print ('total of {} mutations that pass'.format(mut_count), file=sys.stderr)
