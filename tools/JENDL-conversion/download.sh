#!/bin/bash
################################################################################
# Download files and convert Windows newline to Linux newline
################################################################################

wget -O jendlan2005.tar.gz http://wwwndc.jaea.go.jp/ftpnd/ftp/JENDL/jendlan2005.tar.gz
tar -xzf jendlan2005.tar.gz 
mkdir JENDL-AN-2005-linux > /dev/null
for i in $(ls JENDL-AN-2005/*.dat)
do
sed $'s/\r$//' $i > ${i/2005/2005-linux}
done
