from cyvcf2 import VCF
import argparse
import numpy as np
from pyfaidx import Fasta
from collections import defaultdict
from quicksect import IntervalTree
import itertools
import sys

p = argparse.ArgumentParser()
p.add_argument('--vcf')
p.add_argument('--ref')
p.add_argument('--genes')
p.add_argument('--phylop')
p.add_argument('-dnase')
p.add_argument('-exclude')
p.add_argument('-nmer', type=int, default=3)
p.add_argument('-r')
args = p.parse_args()

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


if args.nmer % 2 == 0:
    print ("not a valid nmer count", file=sys.stderr)
    sys.exit()

vcf = VCF(args.vcf)

HOM_ALT = 3
HET = 1
HOM_REF = 0
UNK = 2

reference = Fasta(args.ref)

n_mutation_combos = (4 ** (int(args.nmer - 1))) * 6

samples = ['EAS', 'AMR', 'AFR', 'EUR', 'SAS']

rc = {'G>T':'C>A', 'G>C':'C>G', 'G>A':'C>T', 'T>C':'A>G', 'T>G':'A>C', 'T>A':'A>T'}
rc_nuc = {'A':'T', 'C':'G', 'T':'A', 'G':'C'}

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

expanded_muts = count_mutations(args, rc, rc_nuc)

mut2idx = dict(zip(expanded_muts, range(len(expanded_muts))))
idx2mut = {v:k for k,v in mut2idx.items()}# = dict(zip(expanded_muts, range(len(expanded_muts))))
mut_count = 0

genic_mutation_counts = {s:np.zeros(n_mutation_combos) for s in samples}
ncdng_mutation_counts = {s:np.zeros(n_mutation_combos) for s in samples}

gene_tree = read_exclude(args.genes)
exclude = None
if args.exclude:
    exclude = read_exclude(args.exclude)
dnase = None
if args.dnase:
    dnase = read_exclude(args.dnase)

phylop=read_exclude(args.phylop)

if args.r:
    vcf = vcf(args.r)

for i,v in enumerate(vcf):
    if i % 50000 == 0: 
        d_frac = np.where(genic_mutation_counts['EUR'] > 0)[0].shape[0] / float(n_mutation_combos)
        print ('done with {} variants, currently on chrom {} (catalogued {} of all {}-mers)'.format(i, v.CHROM, d_frac, args.nmer), file=sys.stderr)

    if exclude is not None and len(exclude[v.CHROM].search(v.start, v.end)) > 0: continue
    
    if v.FILTER not in (None, "PASS"): continue
    if len(v.REF) > 1: continue
    if len(v.ALT) > 1: continue
    if len(v.ALT[0]) > 1: continue
    if v.CHROM in ("chrX", "chrY", "X", "Y"): continue
    if v.CHROM not in [str(n) for n in range(1,23)]: continue

    in_gene = False  
    if len(gene_tree['chr' + v.CHROM].search(v.start, v.end)) > 0:
        in_gene = True

    if args.dnase and len(dnase['chr' + v.CHROM].search(v.start, v.end)) == 0: continue
    if len(phylop['chr' + v.CHROM].search(v.start, v.end)) > 0: continue

    #if not all([float(v.INFO.get('{}_AF'.format(s))) > 0 for s in ['AFR', 'EUR', 'EAS', 'SAS', 'AMR']]): continue
    if any([float(v.INFO.get('{}_AF'.format(s))) > 0.98 for s in ['AFR', 'EUR', 'EAS', 'SAS', 'AMR']]): continue

    n_nucs = None
    if args.nmer == 1:
        n_nucs = 0
    else: n_nucs = int((args.nmer - 1) / 2)

    if args.nmer == 1:
        upstream, dnstream = None, None
    else:
        upstream = str(reference[v.CHROM][v.start - n_nucs:v.start]).upper()
        dnstream = str(reference[v.CHROM][v.end:v.end + n_nucs]).upper()

    assert str(reference[v.CHROM][v.start]).upper() == v.REF

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

    for samp in samples:
        allele_count = int(v.INFO.get(samp + '_AF') * v.INFO.get('AN'))
        if allele_count == 0: continue
        if in_gene:
            genic_mutation_counts[samp][kmer_idx] += 1
        elif not in_gene:
            ncdng_mutation_counts[samp][kmer_idx] += 1

mlist = ','.join(mut2idx.keys())

fh1, fh2, = '1kg.k{}.genic.txt'.format(args.nmer), '1kg.k{}.noncoding.txt'.format(args.nmer)

fh1, fh2 = open(fh1, mode='a'), open(fh2, mode='a')

for f in (fh1, fh2):
    print (','.join(['samp', mlist]), file=f)
for samp in genic_mutation_counts:
    g_vals = list(map(int, list(genic_mutation_counts[samp])))
    g_vals = ','.join(list(map(str, g_vals)))
    n_vals = list(map(int, list(ncdng_mutation_counts[samp])))
    n_vals = ','.join(list(map(str, n_vals)))
    print (','.join([samp, g_vals]), file=fh1)
    print (','.join([samp, n_vals]), file=fh2)
