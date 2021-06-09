#! /usr/bin/env python

import pandas
import numpy as np
import allel
import sys

df_samples = pandas.read_csv(sys.argv[2], sep='\t', index_col='Index')

pop1 = sys.argv[3]
pop2 = sys.argv[4]

subpops = {
    pop1: df_samples[df_samples.Site == pop1].index,
    pop2: df_samples[df_samples.Site == pop2].index,
}

genes = np.genfromtxt('genes.txt', dtype='str')

print(sys.argv[3], '-', sys.argv[4])
print('Gene', 'Fst', 'Numerator', 'Denominator')

for i in genes:
    callset = allel.read_vcf(sys.argv[1], region = i)
    pos = callset['variants/POS']
    gt = allel.GenotypeChunkedArray(callset['calldata/GT'])
    acs = gt.count_alleles_subpops(subpops)
    ac1 = allel.AlleleCountsArray(acs[pop1])
    ac2 = allel.AlleleCountsArray(acs[pop2])
    num, den = allel.hudson_fst(ac1, ac2)
    num_sum = np.sum(num)
    den_sum = np.sum(den)
    if den_sum != 0:
        gene_fst_hudson = num_sum / den_sum
    else:
        gene_fst_hudson = 0
    print(i, gene_fst_hudson, num_sum, den_sum)

