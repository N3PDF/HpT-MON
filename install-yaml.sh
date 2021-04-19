#!/bin/bash

branch=master

YAML=yaml-cpp-${branch}

if [[ -d $YAML ]]
then
    echo "Removing existing directory!!"
    rm -rf $YAML
fi

echo "Download package..."
wget https://github.com/jbeder/yaml-cpp/archive/refs/heads/${branch}.zip >& /dev/null
unzip ${branch}.zip >& /dev/null
rm -r ${branch}.zip

cd $YAML

echo "Configure package..."
mkdir build
cd build
cmake .. >& yaml-install.log
if [[ $? != 0 ]]
then
    echo "Error on compilation, check yaml-install.log!!"
    exit -1
fi

echo "Build package..."
make >& yaml-install.log
if [[ $? != 0 ]]
then
    echo "Error on build, check yaml-install.log!!"
    exit -1
fi

export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$(pwd)
echo "Added $PWD to PKG-CONFIG path."
cd ../../
