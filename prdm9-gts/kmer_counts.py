import pysam
import sys
from collections import defaultdict
import argparse
from peddy import Ped
from itertools import permutations
import numpy as np
import re
import csv

p = argparse.ArgumentParser()
p.add_argument('--bams', nargs='*')
p.add_argument('--motifs')
p.add_argument('--alleles')
p.add_argument('-clipping', action='store_true')
p.add_argument('-o', default='out')
p.add_argument('-k', default=15, type=int)
args = p.parse_args()

def count_kmer(seq, k=2):
    """
    return a list of all kmers in a given sequence
    using a simple sliding window of size k

    >>> count_kmer('AAATACCGC', k=2)
    ['AA', 'AA', 'AT', 'TA', 'AC', 'CC', 'CG', 'GC']
    >>> count_kmer('AAATACCGC', k=3)
    ['AAA', 'AAT', 'ATA', 'TAC', 'ACC', 'CCG', 'CGC']
    """
    mers = []
    idx = 0
    while idx + k <= len(seq):
        mers.append(seq[idx:idx+k])
        idx += 1
    return mers

def compare_kmers(uniq_kmers, read_seq, read_bq, read_name, k=30, min_bq=20):
    """
    compare the kmers in a read with kmers
    unique to each ZnF motif
    """

    # store the motifs with read support
    motifs_with_support = []

    # canonical motifs present in A PRDM9 alleles 
    a_motifs = ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J')

    # we'll store info about reads that don't
    # support an A PRDM9 allele for output later
    non_a_reads = []

    read_kmers = count_kmer(read_seq, k=k)
    read_kmers = set(read_kmers)
    for motif in uniq_kmers:
        motif_kmers = uniq_kmers[motif]
        match_kmers_in_read = set(read_kmers).intersection(motif_kmers)
        if len(match_kmers_in_read) == 0: continue
        # make sure all kmer matches are hi-quality, using
        # BQ of all bases in the matching portion of the read
        all_kmer_hi_qual = True 
        for kmer in match_kmers_in_read:
            re_iter = re.finditer(kmer, read_seq)
            for re_match in re_iter:
                s, e = re_match.start(), re_match.end()
                qual_at_match = list(read_bq[s:e])
                if not all([q >= min_bq for q in qual_at_match]): all_kmer_hi_qual = False
        if not all_kmer_hi_qual: continue
        if motif not in a_motifs:
            non_a_reads.append((read_name, 
                                list(match_kmers_in_read)[0]))
        motifs_with_support.append(motif)

    return motifs_with_support, non_a_reads

def get_kmers(bam, args, uniq_motifs, include_sclip=False):
    """
    iterate over all reads in a provided bam file
    and figure out which PRDM9 motifs have support
    """
    # simple reverse complement dict
    rc = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

    # read in alignment file and only consider
    # region amplified in Berg et al. 2010
    try: alnfile = pysam.AlignmentFile(bam, "rc")
    except ValueError: return [], []
    #iterator = alnfile.fetch("5", 23526096, 23527995)
    iterator = alnfile.fetch("chr5", 23526096, 23527995)

    motifs_with_support = defaultdict(int)
    non_a_reads = []

    for read in iterator:
        # filter on MQ
        if read.mapping_quality < 20: continue
        # get sequence including/excluding soft clipped bases.
        # BQ array also depends on inclusion of clipped bases
        seq = None
        if include_sclip:
            seq = read.query_sequence
            quals = list(read.query_qualities) 
        else:
            seq = read.query_alignment_sequence
            quals = list(read.query_alignment_qualities) 
        if seq is None or 'N' in seq: continue 
        
        mws, nar = compare_kmers(uniq_motifs, seq, quals, read.query_name)

        for m in mws: motifs_with_support[m] += 1
        non_a_reads.extend(nar) 

    return motifs_with_support, non_a_reads 


def get_uniq(motifs):
    """
    for each allele, return a list of the kmer motifs
    that are unique to that allele

    >>> get_uniq({'A':['AAA', 'AAT', 'ATC', 'TCG'], 'B':['AAT', 'ATC', 'TCG', 'ATA']})
    defaultdict(<class 'list'>, {'A': ['AAA'], 'B': ['ATA']})
    """
    uniq_motifs = defaultdict(list)
    for motif in motifs:
        cur_kmers = set(motifs[motif])
        rest_kmers = []
        for m in motifs:
            if m == motif: continue
            rest_kmers.extend(motifs[m])

        rest_kmers = set(rest_kmers)
        uniq_motifs[motif] = list(cur_kmers.difference(rest_kmers))
    return uniq_motifs

def generate_prdm9_refseq(alleles, motifs):
    """
    generate a set of references sequences for all prdm9 alleles

    >>> generate_prdm9_refseq({'A':'ABCC', 'B':'ABCD'}, {'A':'ATC', 'B':'CCC', 'C':'ACG', 'D':'GTA'})
    defaultdict(<class 'str'>, {'A': 'ATCCCCACGACG', 'B': 'ATCCCCACGGTA'})
    """

    allele_structures = defaultdict(str)

    for allele in alleles:
        for motif in list(alleles[allele]):
            allele_structures[allele] += motifs[motif]
    return allele_structures


def main(args):
    a_motifs = ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J')

    # generate a list of all kmers present
    # in each of the PRDM9 ZnF motifs
    raw_motifs = defaultdict(str)
    motifs = defaultdict(list)
    motif_junctions = defaultdict(list)
    with open(args.motifs, 'r') as f:
        fh1 = csv.reader(f, delimiter='\t')
        for l1 in fh1:
            raw_motifs[l1[0]] = l1[1]
            motifs[l1[0]] = count_kmer(l1[1], k=args.k)
            with open(args.motifs, 'r') as f2:
                fh2 = csv.reader(f2, delimiter='\t')
                for l2 in fh2:
                    #if l1[0] + l2[0] in motif_junctions: continue
                    #if l2[0] + l1[0] in motif_junctions: continue
                    motif_junctions[l1[0] + l2[0]] = l1[1] + l2[1]

    raw_alleles = defaultdict(str)
    alleles = defaultdict(list)
    with open(args.alleles, 'r') as f:
        fh = csv.reader(f, delimiter='\t')
        for l in fh:
            raw_alleles[l[0]] = l[1]
            mers = count_kmer(l[1], k=1)
            for perm_length in (2,):
                fixed_perms = []
                perms = permutations(mers, perm_length)
                for p in perms:
                    if all([pp == p[0] for pp in p]): continue
                    fixed_perms.append(p)
                alleles[l[0]].extend(fixed_perms)
    
    import glob
    if 'prdm9.refseq.fa' not in glob.glob('*.fa'):
        refs = open('prdm9.refseq.fa', 'a') 
        prdm9_refseqs = generate_prdm9_refseq(raw_alleles, raw_motifs)
        for allele in prdm9_refseqs:
            print ('>{}\n{}'.format(allele, prdm9_refseqs[allele]), file=refs)

    
    uniq_alleles = get_uniq(alleles) 
    uniq_motifs = get_uniq(motifs)
    uniq_motif_junctions = get_uniq(motif_junctions)

    #print ([a for a in uniq_alleles['L17'] if a in uniq_alleles['L21']])
    for i,bam in enumerate(args.bams):
        ev, alt_reads = get_kmers(bam, args, uniq_motifs, include_sclip=args.clipping)
        if ev is None: continue
        sample = bam.split('/')[-1].split('.')[0]

        alleles_with_support = []
        for allele in uniq_alleles:
            alleles = uniq_alleles[allele]
            if len(alleles) == 0: continue
            for a in alleles:
                if all([x in ev for x in a]): alleles_with_support.append(allele)
        
        if len(alleles_with_support) > 0: 
            if len(set(alleles_with_support)) > 1: prdm9 = "unk"
            #else: prdm9 = alleles_with_support[0]
            else: prdm9 = 'non-A'#alleles_with_support[0]

        elif all([k in a_motifs for k in ev]): prdm9 = 'A'

        else:
            non_a_cts = [ev[k] for k in ev if k not in a_motifs]
            a_cts = [ev[k] for k in ev if k in a_motifs]
            #if max(non_a_cts) >= min(a_cts): prdm9 = 'non-A'
            if max(non_a_cts) >= 1: prdm9 = 'non-A'
            #if max(non_a_cts) >= 5: prdm9 = 'non-A'
            else: prdm9 = 'unk'
        #print (ev) 
        print ('\t'.join([sample, 
                            prdm9, 
                            ','.join([k + ':' + str(ev[k]) for k in ev])]))
                            #','.join([r[0] + ':' + r[1] for r in alt_reads])]))

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    main(args)
