import sys
import csv
from cyvcf2 import VCF

vcf = VCF(sys.argv[2])

outfile = open('mouse.cc.private.vcf', 'a')

print (vcf.raw_header.rstrip(), file=outfile)

with open(sys.argv[1], 'r') as fh:
    f = csv.reader(fh)
    for i,l in enumerate(f):
        if i == 0: continue
        for v in vcf(l[1] + ':' + str(int(l[2]) + 0) + '-' + str(int(l[2]) + 0)):
            print (str(v).rstrip(), file=outfile)

