#!/usr/bin/python

import sys
from cyvcf2 import VCF

dat = dict()

# Parse the VCF file to get the list of chromosomes, as well as the position of the first and last SNPs for each
for record in VCF(sys.argv[1]):
  chrom = record.CHROM
  posit = record.POS
  if chrom in dat.keys():
    dat[chrom][1] = posit
  else:
    dat[chrom] = [posit, posit]

# Now generates the tab file:
for chrom, positions in dat.items():
  print("%s\t%i\t%i\t%i\t%i\t%i\t%i" % (chrom, positions[0], positions[1], 0, positions[1] - positions[0], positions[0], positions[1]))

