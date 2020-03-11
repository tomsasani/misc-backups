import sys
import toolshed as ts
import csv

header = "##fileformat=VCFv4.1\n"
header += '\t'.join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"])

with open(sys.argv[1], 'r') as f:
    fh = csv.reader(f, delimiter='\t')
    for i,l in enumerate(fh):
        if i == 0:
            samples = l[4:]
            header += '\t'
            header += '\t'.join(samples)
            print (header)
        else:
            gts = list(map(float, l[4:]))
            gts_f = [g for g in gts if g in (0,1)]
            af = sum(gts_f) / float(len(gts_f))
            ac = sum(gts_f)
            an = len(gts_f)
            chrom = l[0]
            pos = l[1]
            ident = '.'
            ref, alt = l[2], l[3]
            qual = '1'
            filt = 'PASS'
            info = 'AF={};AC={};AN={}'.format(af, int(ac), int(an))
            form = 'GT'
            gt_array = ['/'.join([str(int(g)), str(int(g))]) if g in (0,1) else './.' for g in gts]
            print ('\t'.join([chrom,
                              pos,
                              ident,
                              ref,
                              alt,
                              qual,
                              filt,
                              info,
                              form,
                              '\t'.join(gt_array)]))






