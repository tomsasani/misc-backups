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

mut_count = 0
for i,v in enumerate(vcf):
    if i % 10000 == 0: print ('done with {} variants, on chrom {}'.format(i, v.CHROM), file=sys.stderr)
    if v.FILTER not in (None, "PASS"): continue

    if len(v.REF) > 1: continue
    if len(v.ALT) > 1: continue
    if len(v.ALT[0]) > 1: continue
    if v.INFO.get('MQ') is not None and int(v.INFO.get('MQ')) < 30: continue
    if v.CHROM == "X": continue

    if v.INFO.get('ANN') is None: continue

    gts = v.gt_types

    call_rate = np.sum(gts != 2) / float(len(smp2idx))

    if call_rate < 0.99: continue

    annotation = v.INFO.get('ANN')

    if 'intron' in annotation: continue
    if not 'missense' in annotation: continue
    gene = annotation.split('|')[3]

    n2_lsj1_idx = np.array([v for k,v in smp2idx.items() if k in ('N2', 'LSJ1')])

    het_gts = np.where(gts == HOM_ALT)[0]

    if len(het_gts) > 2: continue

    shared_elems = np.any(np.in1d(het_gts, n2_lsj1_idx))

    if bool(shared_elems) is False: continue

    print ('\t'.join([str(v.CHROM), str(v.start), str(v.end), ','.join([str(x) for x in n2_lsj1_idx]), annotation]))
