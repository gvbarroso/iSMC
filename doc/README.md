# A starter guide for creating the input files for iSMC

`ismc` basic functionally is built around the bio++ libraries (https://github.com/BioPP)
As such, input parameters are specified in an 'options file'. Let's call this options file opt.bpp. To run 'ismc' from the command-line:

```
ismc params=opt.bpp
```

We now describe the options file for `ismc v1`. Although under the hood there are a number of options that enable flexibility for specifying different models and ways to filter the input data, these have default parameters that we don't have to worry about here. The options that we want to specify in opt.bpp are the following (specify them in your preferred order):

## Input/Output file options

These are straightforward options related to the user interface, and NOT to the actual population genetics model.

### dataset label

A STRING parameter that will be appended to the output files, with the purpose of aid organization. For example:

```
dataset_label = altai_neanderthal
```

### input file type

The type of the input sequence data (STRING). It can be either "FASTA", "VCF" or "gVCF".
If "FASTA", 'N' characters and gaps will be masked out by `ismc`. 
If "VCF", an accompanying mask file (in FASTA or BED format, see below) should be provided (see below).
If "gVCF", absent positions will be masked by `ismc`.

```
input_file_type = # FASTA, VCF or gVCF
```

### sequence compression type

The type of compression employed in the sequence input file (STRING). It can be "none" (DEFAULT), "gzip" or "bgzip"

```
seq_compression_type = # none, gzip or bgzip
```

### sequence file path

The relative path to the input sequence file (STRING). For example:

```
seq_file_path = ../data/my_data.vcf.gz
```

### mask file type

The format of the mask file, mandatory if the input file type is VCF (STRING). It can be "FASTA" or "BED". The mask must contain only sites present in the sequence file, and in the exact same order. In the "BED" case, `ismc` will mask out (convert to missing data) sites that are present in this file (i.e., a "negative" mask). In the "FASTA" case, `ismc` will keep mask out sites that are NOT represented by either '1' or 'P' (the callable code).

```
mask_file_type = # FASTA or BED
```

### mask compression type

The type of compression employed in the input file (STRING). It can be "none" (DEFAULT), "gzip" or "bgzip".

```
mask_compression_type = # none, gzip or bgzip
```

### mask file path

The relative path to the input mask file (STRING). For example:

```
mask_file_path = ../data/my_mask.bed
```

### tab file

The relative path to a tab-separated (.tsv) file specifying one genomic 'block' per line (e.g. chromosome or scaffold). Such blocks must mirror the segment structure of the input sequence file.
This file has 7 columns, that must represent the following information:

1: block ID (eg chr1)
2: the start coordinate of the block mapped to the reference genome
3: the end coordinate of the block mapped to the reference genome

NOTE: start and end coordinates are relative to the block, meaning it restarts e.g. at every chromosome.
 
4: '0' for all blocks
5: the difference betwee 3rd & 2nd columns

6: bottom cut-off position from where `ismc_mapper` should start binning single-nucleotide rates into larger genomic windows.
7: top cut-off position until where `ismc_mapper` should start binning single-nucleotide rates into larger genomic windows.

The last two columns are convenient to 'synchronise' coordinates when using `ismc_mapper`. This can happen when you want to consider a range of sites that matches that of another file. For example, if your input sequence data for chromosome 1 starts at position 5,000 and goes until position 10,000,000, and you want to 'synchronise' it with an experimental genetic map that ranges from position 10,000 to 8,000,000, the 6th and 7th columns of the first line of your tab_file should be  9,999 and 7,999,999. Otherwise, if you want to include all sites, they should be the same as the 2nd and 3rd columns.

An example file is given in the 'example' directory.
Specify it in the options file with e.g.:

```
tab_file_path = my_tab.tsv
```

### order of haplotype indices for building pseudo-diploids

When data is phased and multiple haplotypes are available, iSMC can combine them in user-defined pairs. These are specified as comma-separated INTEGERS enclosed by parentheses, indexed from zero (DEFAULT = (0,1) implies that only a t).
For example, to arrange three haplotypes into two pairs of genomes, where the first pair is made up of haplotypes #0 and 1 and the second pair is made up of haplotypes #1 and #2:

```
diploid_indices = (0,1,1,2) // DEFAULT = (0,1)
```

number_threads = 8
time_steps = $(T)
factor_order = $(O)
demes_file = $(F)
```

we can execute e.g.:

```
momentspp params=opt.bpp T=2000 O=75 F=model_pop1.yaml
```

Once again, an example of this functionality is given in two_locus_time.Rmd
