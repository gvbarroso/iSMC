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
7: top cut-off position until where `ismc_mapper` should bin single-nucleotide rates into larger genomic windows.

The last two columns are convenient to 'synchronise' coordinates when using `ismc_mapper`. This can happen when you want to consider a range of sites that matches that of another file. For example, if your input sequence data for chromosome 1 starts at position 5,000 and goes until position 10,000,000, and you want to 'synchronise' it with an experimental genetic map that ranges from position 10,000 to 8,000,000, the 6th and 7th columns of the first line of your tab_file should be  9,999 and 7,999,999. Otherwise, if you want to include all sites, they should be the same as the 2nd and 3rd columns.

An example file is given in the 'example' directory.
Specify it in the options file with e.g.:

```
tab_file_path = my_tab.tsv
```

### order of haplotype indices for building pseudo-diploids

When data is phased and multiple haplotypes are available, iSMC can combine them in user-defined pairs. These are specified as comma-separated INTEGERS enclosed by parentheses, indexed from zero. iSMC will then combine these indices in non-overlapping pairs.
For example, to arrange three haplotypes into two pairs of genomes, where the first pair is made up of haplotypes #0 and 1 and the second pair is made up of haplotypes #1 and #2:

NOTE: DEFAULT=?

```
diploid_indices = (0,1,1,2) // DEFAULT = ?
```

### number of computing threads

The number of computing threads `ismc` is allowed to use (INTEGER), with DEFAULT equal to all available threads in the machine. 
Note that the typical axes of parallelization possible in `ismc` are over haplotype pairs ("diploids") and over genomic blocks (usually chromosomes).
Also note that parallelization in the decoding step often incurs in large memory use (see below).

```
number_threads = 4 # DEFAULT = all threads
```

### iSMC task 1

Should `ismc` perform parameter optimization (BOOLEAN)? 

```
optimize = true
```

### iSMC task 2

Should `ismc` perform posterior decoding (BOOLEAN)?

```
decode = false
```

Note: As the optimization and decoding steps require different computational resources, it is often a good idea to run `ismc` separately for each of these tasks.

### Recovering an interrupted optimization

Depending on the model complexity, sample size and genome length, the optimization step can take a long time to finish. This means that there are a lot of opportunities for something to go wrong in the middle of optimization (e.g., you reach the run-time limit you requested in the computing cluster). For this reason, `ismc` keeps an updated file with intermediate values for best fit parameters, called 'backup_params.txt'. You can then resume an interrupted optimization with the following option (BOOLEAN).

```
resume_optim = true # DEFAULT = false
```

### Speeding up posterior deconding 1

Should `ismc` perform posterior decoding in parallel (over pairs of haplotypes, BOOLEAN)?

```
decode_diploids_parallel = true # DEFAULT = false
```

### Speeding up posterior deconding 1

Should `ismc` perform posterior decoding in parallel (over genomic blocks, BOOLEAN)?

```
decode_breakpoints_parallel = true # DEFAULT = false
```

Note that memory use increases linearly with paralellizing over diploids and/or blocks.

### Decoding fragment size

`ismc` should not decode entire chromosomes at once as this would consume too much memory and the process would most certainly be killed. Therefore, we can specify the size of fragments to decode (INTEGER), effectively breaking up the task (and resulting in several output files that will then be collected by `ismc_mapper`).

```
fragment_size = 1000000 # DEFAULT = 3000000
```

## Coalescent model options

The following options are more intimately related to the specification of the evolutionary model.

### Number of time intervals

The number of discretized time intervals in the TMRCA distribution (INTEGER). This is related to the MSMC2 option, except `ismc` uses a splines model to reduce the number of model parameters. Up to a point, more time intervals means more power to detect genealogical transitions, but also a larger number of hidden states in the HMM, implying longer run-time and memory use. In practice, we have found 30 to often represent a decent trade-off.

```
number_intervals = 30 # DEFAULT = 40 

```

### Number of splines knots

The number of "knots" where to anchor the splines curves taming the coalescence rates over time (INTEGER). These knots are distributed in the temporal axis, and define the number of curves to fit (which os equal to the number of knots + 1). By default, they are constrained to be strictly positive, to have derivative equal to zero at each end and to meet the adjacent splints at the junction points (Inverse coalescence rates are then mapped from the splines curves at each defined time interval).

If this writing seems to convoluted, note that in practice having between 2 and 3 knots is usually a good choice.

In principle it is possible to specify an initial and a maximum number of knots, in which case `ismc` will fit progressively more complex models and select among them, but this is probably overkill. That is, simply specify both init and max numbers to be the same.

```
init_number_knots = 2 # DEFAULT = 3
max_number_knots = 2 # DEFAULT = 3

```

### Number of discretized recombination categories

The number of categories into which to discretize the (Gamma) distribution of recombination rates (INTEGER). The default value of 1 means no recombination rate variation along the genome.

```
number_rho_categories = 5 # DEFAULT = 1
```

### Number of discretized mutation categories

The number of categories into which to discretize the (Gamma) distribution of mutation rates (INTEGER). The default value of 1 means no mutation rate variation along the genome.

```
number_theta_categories = 5 # DEFAULT = 1
```

## Putting it all together

The options file supports passing bash-style variables as input. For example, if the options file (opt.bpp) looks like the following

```
dataset_label = $(DATA)_gamma5_$(NB_KNOTS)knots

sequence_file_path = $(DATA).vcf.gz
input_file_type = VCF
mask_file_path =  $(DATA)_negmask.bed.gz 
mask_file_type = BED
seq_compression_type = gzip
mask_compression_type = gzip
tab_file_path = $(DATA).tab
diploid_indices = (0,1)

init_number_knots = $(NB_KNOTS)
max_number_knots = $(NB_KNOTS)

optimize = true 
resume_optim = false

decode = false
fragment_size = 10000000

number_theta_categories = 1
number_rho_categories = 5
number_ne_categories = 1
number_intervals = 30

```

we can execute, e.g.:

```
ismc params=opt.bpp DATA=altai_chr1 NB_KNOTS=3
```

