#!/bin/bash

VERSION=6.2.3

rm -r lhapdf6 >& /dev/null
mkdir -p lhapdf6
cd lhapdf6
echo "Installing LHAPDF version $VERSION"

rm pdfinstall.log >& /dev/null
rm -r bin include lib share LHAPDF-${VERSION} >& /dev/null

echo "Downloading tar file..."
wget http://www.hepforge.org/archive/lhapdf/LHAPDF-${VERSION}.tar.gz
tar -xzvf LHAPDF-${VERSION}.tar.gz
rm LHAPDF-${VERSION}.tar.gz

cd LHAPDF-${VERSION}
echo "Configuring..."
export CPPFLAGS="-P"
./configure --prefix=`cd .. && pwd` --disable-python

if [[ $? != 0 ]]
then
    echo "Error on configuration, check pdfinstall.log"
    exit -1
fi

echo "Compiling..."
make >& pdfinstall.log
if [[ $? != 0 ]]
then
    echo "Error on compilation, check pdfinstall.log"
    exit -1
fi

echo "Installing..."
make install >& pdfintall.log
if [[ $? != 0 ]]
then
    echo "Error on installation, check pdfinstall.log"
    exit -1
fi
cd ..

echo "LHAPDF software installation completed"

echo "Downloading PDF sets..."
pdflist=`cat ../pdfsets.list`

cd share/LHAPDF
for pdf in $pdflist
do
    pdf=`echo $pdf | cut -d. -f1`
    if [[ -e $pdf ]]
    then
	continue
    fi
    wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/${pdf}.tar.gz >& /dev/null
    if [[ $? != 0 ]]
    then
	echo "$pdf download failed"
	continue
    fi
    tar -xzvf ${pdf}.tar.gz >& /dev/null
    echo "Extraction complete!"
    rm ${pdf}.tar.gz
done

cd ../..
cd LHAPDF-${VERSION}
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$(pwd)
echo "Added $PWD to PKG-CONFIG path."
cd ../../
