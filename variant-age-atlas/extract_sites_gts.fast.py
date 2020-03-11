from cyvcf2 import VCF
import argparse
import numpy as np
from pyfaidx import Fasta
from collections import defaultdict
from quicksect import IntervalTree
import itertools
import sys
import time
import csv
import gzip

p = argparse.ArgumentParser()
p.add_argument('--csv')
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
reference = Fasta(args.ref)

# count numbers of possible mutations
n_mutation_combos = (4 ** (int(args.nmer - 1))) * 6

rc = {'G>T':'C>A', 'G>C':'C>G', 'G>A':'C>T', 'T>C':'A>G', 'T>G':'A>C', 'T>A':'A>T'}
rc_nuc = {'A':'T', 'C':'G', 'T':'A', 'G':'C'}

expanded_muts = count_mutations(args, rc, rc_nuc)

mut2idx = dict(zip(expanded_muts, range(len(expanded_muts))))
idx2mut = {v:k for k,v in mut2idx.items()}
mut_count = 0

age_array = np.arange(0,200000,5000)

ncdng_mutation_counts = np.zeros((n_mutation_combos, age_array.shape[0]))

r_count = 0
with (gzip.open(args.csv, 'rt') if args.csv.endswith('gz') else open(args.csv, 'r')) as fh:
    f = csv.reader(fh)
    for i,l in enumerate(f):
        l = [x.rstrip()[1:] for x in l]
        if 'ID' in l[0]: continue
        if i % 1000 == 0: print ('done with {} variants'.format(i))
       
        chrom = l[1]
        start = int(l[2]) - 1
        end = int(l[2])

        n_nucs = None
        if args.nmer == 1:
            n_nucs = 0
        else: n_nucs = int((args.nmer - 1) / 2)

        if args.nmer == 1:
            upstream, dnstream = None, None
        else:
            upstream = str(reference[chrom][start - n_nucs:start]).upper()
            dnstream = str(reference[chrom][end:end + n_nucs]).upper()

        ref, alt = l[3], l[4]

        assert str(reference[chrom][start]).upper() == ref

        if upstream is not None and dnstream is not None:
            if not (len(str(upstream)) == n_nucs and len(str(dnstream)) == n_nucs): continue


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

        v_age = float(l[11])
        v_age_idx = np.histogram(v_age, bins=age_array)
        v_age_idx = np.where(v_age_idx[0] > 0)[0]

        ncdng_mutation_counts[kmer_idx][v_age_idx] += 1

header_cols = ','.join(['-'.join([str(age_array[n-1]), str(age_array[n])]) for n,_ in enumerate(age_array) if n > 0])
header_rows = [*mut2idx]

fh = 'k{}.{}.age.mutation_spectra.txt'.format(args.nmer, args.outname)

fh = open(fh, mode='a')
print (','.join(['mutation_type', header_cols]), file=fh)
for i,samp in enumerate(ncdng_mutation_counts):
    n_vals = list(map(int, list(samp)))
    n_vals.insert(0, header_rows[i])
    n_vals = ','.join(list(map(str, n_vals)))
    print (n_vals, file=fh)
