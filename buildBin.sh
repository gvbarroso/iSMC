#! /bin/sh
arch=`uname -m`
version=1.0.0-1

strip src/ismc
strip src/ismc_mapper
strip src/ismc_mapper2
tar cvzf ismc-${arch}-bin-static-${version}.tar.gz src/ismc src/ismc_mapper src/ismc_mapper2

