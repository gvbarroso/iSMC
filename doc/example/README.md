# Test example

Data from one individual from the YRI population (1000 Genomes), chromosome 9 only.

## Data preparation

Generate a default tab file:

```bash
python3 ../../tools/Python/generate_tab.py chr9.YRI.NA18486.norm.vcf.gz
```

Extract mask file for chr9:

```bash
zgrep chr9 20160622.allChr.mask.bed.gz | gzip -c > 20160622.chr9.mask.bed.gz 
```bash

Create a negative mask:

```bash
echo "chr9	1	138394717" > chr9.bed   
bedtools subtract -a chr9.bed -b 20160622.chr9.mask.bed.gz | gzip -c > 20160622.chr9.negmask.bed.gz
```


## Estimate parameters

We compare spline models with 2, 3 or 4 knots:

```bash
ismc param=opt.bpp
```

The model with four knots is selected.

## Run the posterior decoding


```bash
ismc param=opt.bpp resume_optim=true optimize=false decode=true
```

Then we try with a different discretization:

```bash
ismc param=opt.bpp resume_optim=true optimize=false decode=true dataset_label=example2 "rho_boundaries=(0,0.0001,0.001,0.01,0.1,1,10)"
```

Data retrived from https://www.simonsfoundation.org/simons-genome-diversity-project/ on Aug 10 2018

Original VCF file: LP6005592-DNA_D03.annotated.nh2.variants.vcf.gz
Original mask file: x75.fa.gz

Both files were filtered and only chromosomes 1-4 were kept.
To follow through this example run, execute the following steps:

1. Install the programs ismc and ismc_mapper, or copy executable files to this directory.

2. Run ismc from the command line using the params file provided:
    ismc params=opt.bpp

   NOTE: you can modify the opt.bpp file in order to better suit your computational resources.

3. ismc will output several files. To obtain binned recombination maps, run ismc_mapper from the command line using the params file provided:
    ismc_mapper params=map.bpp

