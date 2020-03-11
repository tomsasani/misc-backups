import csv
import sys
from peddy import Ped
from collections import defaultdict

samples = [s for s in Ped(sys.argv[2]).samples()]
kids = [s for s in samples if s.mom is not None and s.dad is not None]

gt = defaultdict()
alleles = defaultdict()
reads = defaultdict()

with open(sys.argv[1]) as f:
    fh = csv.reader(f, delimiter='\t')
    for l in fh:
        gt[l[0]] = l[1]
        alleles[l[0]] = l[2]
        reads[l[0]] = l[3]

md_seen = []

inh_errors = 0
total_inh = 0

for s in kids:
    if s.dad.sample_id == '8477': continue
    parent_gts = [gt[p] for p in (s.mom.sample_id, s.dad.sample_id)]
    #if 'unk' in parent_gts: continue
    #if gt[s.sample_id] == 'unk': continue
    if gt[s.sample_id] not in parent_gts: 
        print ('\t'.join([s.family_id, s.sample_id, gt[s.sample_id], gt[s.mom.sample_id], gt[s.dad.sample_id], alleles[s.sample_id], alleles[s.mom.sample_id], alleles[s.dad.sample_id], reads[s.sample_id]]))
    sibs = [sa.sample_id for sa in kids if sa.mom.sample_id == s.mom.sample_id]
    sib_gts = [gt[samp] for samp in sibs]

    if s.mom.sample_id + '_' + s.dad.sample_id not in md_seen:

        total_inh += len(sibs)
        inh_errors += len([gt for gt in sib_gts if gt not in parent_gts])
    md_seen.append(s.mom.sample_id + '_' + s.dad.sample_id)

print (inh_errors, total_inh)

