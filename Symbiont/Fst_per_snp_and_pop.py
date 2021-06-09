#! /usr/bin/env python

import pandas
import numpy as np
import allel
import sys

np.seterr(divide='ignore', invalid='ignore')

callset = allel.read_vcf(sys.argv[1])

gt = allel.GenotypeChunkedArray(callset['calldata/GT'])

chrom = callset['variants/CHROM']
pos = callset['variants/POS']

df_samples = pandas.read_csv(sys.argv[2], sep='\t', index_col='Index')

pop1 = sys.argv[3]
pop2 = sys.argv[4]

subpops = {
    pop1: df_samples[df_samples.Site == pop1].index,
    pop2: df_samples[df_samples.Site == pop2].index,
}

acs = gt.count_alleles_subpops(subpops)

ac1_all = allel.AlleleCountsArray(acs[pop1])
ac2_all = allel.AlleleCountsArray(acs[pop2])

num, den = allel.hudson_fst(ac1_all, ac2_all)
snp_fst_hudson = num / den

print(sys.argv[3], '-', sys.argv[4])
print('CHROM', 'POS', 'Fst', 'Numerator', 'Denominator')
for a, b, c, d, e in zip(chrom, pos, snp_fst_hudson, num, den):
    print(a, b, c, d, e)

acu = allel.AlleleCountsArray(acs[pop1][:] + acs[pop2][:])
flt = acu.is_segregating() # & (acu.max_allele() == 1)

ac1 = allel.AlleleCountsArray(acs[pop1].compress(flt, axis=0)[:, :2])
ac2 = allel.AlleleCountsArray(acs[pop2].compress(flt, axis=0)[:, :2])

fst_hudson, se_hudson, vb_hudson, vj_hudson = allel.average_hudson_fst(ac1, ac2, blen=100)

with open("Fst_per_pop.txt", "a") as out:
    print(sys.argv[3], '-', sys.argv[4], ': %.04f +/- %.04f' % (fst_hudson, se_hudson), file=out)


