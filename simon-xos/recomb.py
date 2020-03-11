import os
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1" 
from cyvcf2 import VCF
from peddy import Ped
import numpy as np
import sys
import pickle
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from quicksect import IntervalTree
import doctest
import time
import gzip

p = argparse.ArgumentParser()
p.add_argument('--vcf')
p.add_argument('--ped')
p.add_argument('-exclude')
#p.add_argument('-kid', default='1016')
p.add_argument('-chrom')
p.add_argument('-ab', type=float, default=0.3)
p.add_argument('-inf_parent', default=None)
p.add_argument('-chr_prefix', action="store_true")
args = p.parse_args()

HET, HOM_REF, HOM_ALT = 1, 0, 3
UNKNOWN = -2

def ab_is_good(i, ad, rd, gt=1, het_ab=0.3):
	min_ab = {0:0., 1:het_ab, 3:0.95}
	max_ab = {0:0., 1:1-het_ab, 3:1.0}

	ab = ad[i] / float(rd[i] + ad[i])
	if ab >= min_ab[gt] and ab <= max_ab[gt]: return True
	return False

def read_exclude(exclude, chrom=None):
    """
    generate a binary search tree
    from a provided exclude file
    """
    if exclude is None: return None
    tree = defaultdict(IntervalTree)
    if chrom is not None:
        chrom = chrom.split(':')[0]
    added = 0
    for i, line in enumerate((gzip.open if exclude.endswith(".gz") else open)(exclude)):
        toks = line.rstrip().split('\t')
        if i == 0:
            try: int(toks[1])
            except ValueError: continue
        if chrom is not None and toks[0] != chrom: continue
        added += 1
        tree[toks[0]].add(int(toks[1]), int(toks[2]))
        if added == 0:
            sys.stderr.write("didnt add any intervals to exclude for chrom:%s\n" % chrom)
    return tree

def adjust_unknown(irow, row, kids):
    """
    if a sample at a given sites has a bad genotype
    (i.e., low depth, genotype quality, allele balance),
    assign their genotype at that site to their genotype
    at the previous site

    >>> adjust_unknown(2, np.array([0, 1, 1, -1]), np.array([[0, 1, 1, 0],[1, 0, 0, 1],[0, 1, 1, -1], [0, 1, 1, 1]]))
    array([0, 1, 1, 0])
    >>> adjust_unknown(1, np.array([1, -1, 1, -1]), np.array([[0, 1, 1, 0],[1, -1, 1, -1],[0, 1, 1, 0], [0, 1, 0, 0]]))
    array([1, 1, 1, 0])
    >>> adjust_unknown(1, np.array([1, -1, 0, -1, 0]), np.array([[0, 1, 1, 0, 1], [1, -1, 0, -1, 0]]))
    array([1, 0, 0, 1, 0])
    >>> adjust_unknown(1, np.array([0, -1, 0, -1, 0]), np.array([[0, 1, 1, 0, 1], [0, -1, 0, -1, 0]]))
    array([0, 0, 0, 1, 0])
    """

    unkn = np.where(row == -1)
    kn = np.where(row != -1)
    maj_switch = False
    switch_idx = row[kn] == kids[irow - 1][kn]
    if sum(switch_idx) >= kn[0].shape[0] / 2:#kn[0].shape[0] - 1:
        maj_switch = True
    if maj_switch:
        row[unkn] = kids[irow - 1][unkn]
    else:
        row[unkn] = 1 - kids[irow - 1][unkn]

    return row

def contiguous_regions(condition):
    """Finds contiguous True regions of the boolean array "condition". Returns
    a 2D array where the first column is the start index of the region and the
    second column is the end index."""

    # Find the indicies of changes in "condition"
    d = np.diff(condition)
    idx, = d.nonzero() 

    # We need to start things after the change in "condition". Therefore, 
    # we'll shift the index by 1 to the right.
    idx += 1

    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]

    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size] # Edit

    # Reshape the result into two columns
    idx.shape = (-1,2)
    return idx

def xo(states, inf_sites, chrom=None, parent=None):
    """
    find crossovers in a matrix of states x positions
    """

    xos = defaultdict()
    for parent in states:
        p_states = states[parent]
        p_inf_sites = inf_sites[parent]
        result = np.zeros(p_states.shape, dtype=bool)

        output_result = np.equal(p_states[:,0], p_states[:,1])

        t_blocks = contiguous_regions(output_result)

        threshold = 50 

        for i,b in enumerate(t_blocks):
            s, e = b[0], b[1]

            if e - s > threshold:
                if i == 0:
                    xo_start, xo_end = 0, s
                else:
                    xo_start = t_blocks[i-1][-1]
                    xo_end = s


    return xos

def is_good_site(idx, quals, gts, ad, rd, het_ab=0.3):
    if quals[idx] < 20: return False
    if ad[idx] + rd[idx] < 10: return False
    if gts[idx] == UNKNOWN: return False
    if not ab_is_good(idx, ad, rd, gt=gts[idx], het_ab=het_ab): return False
    return True

def get_xo(args, haplotypes, inf_positions):
    if args.chrom is None:
        if args.chr_prefix is None:
            chroms = map(lambda x: str(x), range(1,22))
            chroms.append('X')
        else:
            chroms = map(lambda x: str(x), range(1,22))
            chroms = ['chr' + c for c in chroms]
            chroms.append('chrX')
        for chrom in chroms:
            xos = xo(haplotypes, inf_positions, sibs, chrom=chrom, parent=args.parent)
            #plot(xos, haplotypes, sibs, inf_positions, args)
        return None
    elif args.chrom:
        xos = xo(haplotypes, inf_positions, chrom=args.chrom)
        plot(xos, haplotypes, inf_positions)
        return xos

def plot(xos, haplotypes, inf_positions):

    for p in xos:
        p_xos = xos[p].transpose()

        cmap = "Reds"# if args.parent == "dad" else "Blues"
        p_inf_positions = np.array(inf_positions[p])
        for kid in p_xos:
            f, ax = plt.subplots(figsize=(10,10))
            color_dict = {False:'grey', True:'firebrick'}
            f_yval = lambda x: color_dict[x]
            yvals = np.vectorize(f_yval)
            colors = yvals(kid)
            state_0 = np.where(colors == 'grey')[0]
            state_1 = np.where(colors == 'firebrick')[0]
            if state_0.shape[0] == 0:
                state_1_y = np.ones(state_1.shape[0])
                ax.scatter(p_inf_positions[state_1], state_1_y, facecolors='none', edgecolors='firebrick', alpha=0.5)
            elif state_1.shape[0] == 0:
                state_0_y = np.zeros(state_0.shape[0])
                ax.scatter(p_inf_positions[state_0], state_0_y, facecolors='none', edgecolors='grey', alpha=0.5)
            else:
                state_1_y = np.ones(state_1.shape[0])
                state_0_y = np.zeros(state_0.shape[0])
                ax.scatter(p_inf_positions[state_1], state_1_y, facecolors='none', edgecolors='firebrick', alpha=0.5)
                ax.scatter(p_inf_positions[state_0], state_0_y, facecolors='none', edgecolors='grey', alpha=0.5)
            
            ax.set_yticks([1])
            #ax[i].set_yticklabels([sibs[i].sample_id])
            ax.tick_params(axis='y', labelrotation=90)
            ax.set_ylim(-0.5, 1.5)
            sns.despine(ax=ax)
            

            ax.set_xlabel('Genomic position (bp)')
            f.savefig(p + '_' + args.chrom + '_' + str(args.ab) + '_' + 'xo.png', dpi=300, bbox_inches="tight")
        
def run(args):
    print ('\t'.join(['chrom', 'start', 'end', 'sample_id', 'parent_id', 'n_vars', 'hap_start']))
    vcf = VCF(args.vcf)
    ped = Ped(args.ped)

    samples = [s for s in ped.samples()]
    kids = [s for s in samples if s.mom is not None and s.dad is not None]
    fams = set([s.family_id for s in samples])

    smp2ped = dict(zip([s.sample_id for s in samples], samples))

    exclude = read_exclude(args.exclude)

    #kid = smp2ped[args.kid].sample_id

    #dad, mom = smp2ped[kid].paternal_id, smp2ped[kid].maternal_id

    #samples_in_ped = [s.sample_id for s in samples]
    #samples_in_fam = [s.sample_id for s in samples if s.family_id == smp2ped[kid].family_id]

    # restrict VCF to the samples in the current family
    #vcf.set_samples(samples_in_ped)

    smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))

    bad_positions = []
    v_feats = []

    haps = defaultdict()
    inf_positions = defaultdict(list)

    fam_dict = defaultdict(lambda: defaultdict())

    
    for fam in fams:
        mom = [s for s in samples if s.family_id == fam and s.mom is None and s.dad is None and s.sex == 'female'][0]
        dad = [s for s in samples if s.family_id == fam and s.mom is None and s.dad is None and s.sex == 'male'][0]
        sibs = [s.sample_id for s in kids if s.dad == dad and s.mom == mom]

        fam_dict[fam]['mom'] = mom
        fam_dict[fam]['dad'] = dad
        fam_dict[fam]['sibs'] = ','.join(sibs) 

    nused, i, report_at, t0 = 0, 0, 1000, time.time()
    for i,v in enumerate(vcf(args.chrom)):
        if i != 0 and i % report_at == 0:
            persec = i / float(time.time() - t0)
            print("%s:%d (%.1f/sec) %.2f%% informative (%d/%d variants)" % (v.CHROM, v.POS,
                  persec, 100.0 * nused/i, nused, i), file=sys.stderr)

        if args.exclude and len(exclude[v.CHROM].search(v.start, v.end)) > 0: continue
        if v.var_type != 'snp': continue
        if v.FILTER not in ('PASS', None): continue
        if v.call_rate < 0.90: continue 
        if len(v.ALT) > 1: continue
        if len(v.ALT[0]) > 1: continue
        
        gts = v.gt_types
        quals = v.gt_quals
        rd, ad = v.gt_ref_depths, v.gt_alt_depths

        for fam in fam_dict:

            mom = fam_dict[fam]['mom']
            dad = fam_dict[fam]['dad']
            sibs = fam_dict[fam]['sibs'].split(',')

            if args.inf_parent and args.inf_parent not in (mom, dad): continue

            try:
                mi, di = smp2idx[mom.sample_id], smp2idx[dad.sample_id]
            except KeyError: sys.exit()

            # ensure we're at an informative site

            if gts[mi] != HOM_REF and gts[di] != HOM_REF: continue
            if gts[mi] == HOM_REF and gts[di] == HOM_REF: continue
            
            sib_gts = [gts[smp2idx[k]] for k in sibs]

            inf_parent = None
            
            if (gts[mi] == HET and gts[di] == HOM_REF): inf_parent = mom.sample_id
            elif (gts[di] == HET and gts[mi] == HOM_REF): inf_parent = dad.sample_id
            else: continue
            
            # check that both parents are "high-quality"
            if not is_good_site(mi, quals, gts, ad, rd, het_ab=args.ab): continue
            if not is_good_site(di, quals, gts, ad, rd, het_ab=args.ab): continue

            if args.inf_parent and inf_parent != args.inf_parent: continue

            # catalog the "states" of each child w/r/t the informative parent
            # if kids are HETs, their state w/r/t to the INF parent is 0
            states = []
            sib_pass = []
            for i,k in enumerate(sibs):
                k_idx = smp2idx[k]
                k_pass = is_good_site(k_idx, quals, gts, ad, rd, het_ab=args.ab)
                sib_pass.append(k_pass)
                k_state = -1
                if v.CHROM in ('chrX', 'X') and k.sex == "male":
                    k_state = 0 if gts[k_idx] == HOM_ALT else 1
                else:
                    k_state = 0 if gts[k_idx] == HET else 1

                states.append(k_state)

            if not all([x is True for x in sib_pass]): continue

            if sum([s < 0 for s in states]) > 0: continue
            
            if inf_parent not in haps:
                haps[inf_parent] = np.array(states)
            else:
                haps[inf_parent] = np.vstack((haps[inf_parent], np.array(states)))

            inf_positions[inf_parent].append(v.start)

        nused += 1

    persec = i / float(time.time() - t0)
    xos = get_xo(args, haps, inf_positions)
#    print (xos)

#if __name__ == "__main__":
doctest.testmod()
run(args)
