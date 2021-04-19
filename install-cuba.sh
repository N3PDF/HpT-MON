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

./configure
if [[ $? != 0 ]]
then
    echo "Error on configuration, check cuba-install.log!!"
    exit -1
fi

make >& cuba-logs.txt
if [[ $? != 0 ]]
then
    echo "Error in build, check cuba-install.log!!"
    exit -1
fi
make install
