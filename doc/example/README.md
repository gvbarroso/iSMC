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

