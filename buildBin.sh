#! /bin/sh
arch=`uname -m`
version=0.0.25-1

strip src/ismc
strip src/ismc_mapper
tar cvzf ismc-${arch}-bin-static-${version}.tar.gz src/ismc src/ismc_mapper

