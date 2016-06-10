#!/bin/bash
################################################################################
# Fix ENDF notation using patches
################################################################################

cp -R JENDL-AN-2005-linux JENDL-AN-2005-linux-patched
cd JENDL-AN-2005-linux-patched
patch -s -p1 <../451-lines.patch
patch -s -p1 <../headline.patch
