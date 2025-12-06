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
```

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

Plot the demography:

```r
demo <- read.table("example_demography.txt", header = TRUE)

require(ggplot2)
p <- ggplot(demo, aes(x=left_time_boundary*20/1.25e-8, y=(1/lambda)/(2*1.25e-8))) +
         geom_step() +
         scale_x_log10() + scale_y_log10()
p
```

## Run the posterior decoding


```bash
ismc param=opt.bpp resume_optim=true optimize=false decode=true
```

We then average estimates in windows:

```bash
ismc_mapper param=map.bpp
```

Then we try with a different discretization (we first copy the parameter estimates):

```bash
cp example_estimates.txt example2_estimates.txt
ismc param=opt.bpp resume_optim=true optimize=false decode=true dataset_label=example2 "rho_boundaries=(0,0.0001,0.001,0.01,0.1,1,10)"
ismc_mapper params=map.bpp dataset_label=example2
```

We then average estimates in windows:

```bash
ismc_mapper param=map.bpp dataset_label=example2
```

Compare the results:

```r
read.data <- function(path, rate, wsize)
{
  dat <- read.table(path, header = TRUE)
  dat$Rate <- rate
  dat$WSize <- wsize
  return(dat)
}
map1.25kb  <- read.data("example.rho.25kb.bedgraph", "Rate 1", "25 kb")
map2.25kb  <- read.data("example2.rho.25kb.bedgraph", "Rate 2", "25 kb")
map1.250kb <- read.data("example.rho.250kb.bedgraph", "Rate 1", "250 kb")
map2.250kb <- read.data("example2.rho.250kb.bedgraph", "Rate 2", "250 kb")
map1.1Mb   <- read.data("example.rho.1Mb.bedgraph", "Rate 1", "1 Mb")
map2.1Mb   <- read.data("example2.rho.1Mb.bedgraph", "Rate 2", "1 Mb")

l <- list()
l[[1]] <- map1.25kb
l[[2]] <- map1.250kb
l[[3]] <- map1.1Mb
l[[4]] <- map2.25kb
l[[5]] <- map2.250kb
l[[6]] <- map2.1Mb
require(data.table)
dat <- rbindlist(l)
dat$WSize <- ordered(dat$WSize, levels = c("25 kb",  "250 kb", "1 Mb"))

p <- ggplot(dat, aes(x=chromStart, y=NA18486)) + geom_line() + facet_grid(Rate~WSize)
p
```

Plot the maps agains each others:

```r
require(tidyverse)
datw <- dat %>% pivot_wider(names_from=Rate, values_from=NA18486)

p <- ggplot(datw, aes(x = .data[["Rate 1"]], y = .data[["Rate 2"]])) + geom_point() + facet_wrap(~WSize) + geom_abline()
p
```

Adding more mass points at the extremes provides better posterior averages.
