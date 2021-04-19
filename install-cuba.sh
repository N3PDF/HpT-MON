#!/bin/bash

version=4.2.1
cuba=Cuba-${version}
cubasrc=`ls $cuba/src`

if [[ x$cubasrc == x ]];
then
    rm -r Cuba-${version}
    wget http://www.feynarts.de/cuba/${cuba}.tar.gz
    tar -xzf ${cuba}.tar.gz
    rm ${cuba}.tar.gz
fi

cd ${cuba}
# sed -i "s/MINSLICE 10/MINSLICE 1/" src/common/Parallel.c
# export CFLAGS="-fPIC -fcommon"
./configure
make
