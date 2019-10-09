#! /bin/sh
arch=`uname -m`
version=0.0.13-1

strip src/ismc
tar cvzf ismc-${arch}-bin-static-${version}.tar.gz --directory=src ismc

