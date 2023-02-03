#!/bin/bash

mkdir tools
echo "Installing tools to run the experiments"
# install sharpSAT
echo "Downloading and compiling sharpSAT"
wget -q --show-progress https://github.com/marcthurley/sharpSAT/archive/v12.08.1.zip 
unzip -q v12.08.1.zip 
mv  sharpSAT-12.08.1 tools/sharpSAT
cd tools/sharpSAT/Release
make
cd ../../

# download conjure
echo "Downloading Conjure"
wget -q --show-progress https://github.com/conjure-cp/conjure/releases/download/v2.3.0/conjure-v2.3.0-linux-solvers.zip
wget -q --show-progress https://github.com/conjure-cp/conjure/releases/download/v2.3.0/conjure-v2.3.0-linux.zip
unzip -q conjure-v2.3.0-linux-solvers.zip
mv conjure-v2.3.0-linux-solvers conjure
unzip -q -j conjure-v2.3.0-linux.zip "conjure-v2.3.0-linux/conjure" -d conjure


# download ASP tools
echo "Downloading ASP tools"
mkdir ASP_tools
cd ASP_tools
wget -q --show-progress https://research.ics.aalto.fi/software/asp/download/binary-x86-64/lp2normal-1.14 
wget -q --show-progress https://research.ics.aalto.fi/software/asp/download/binary-x86-64/lp2sat-1.24
wget -q --show-progress https://research.ics.aalto.fi/software/asp/download/binary-x86-64/lp2normal-2.18
wget -q --show-progress https://research.ics.aalto.fi/software/asp/download/binary-x86-64/lp2atomic-1.17
wget -q --show-progress https://downloads.sourceforge.net/project/potassco/gringo/4.5.4/gringo-4.5.4-linux-x86_64.tar.gz
tar -zxvf gringo-4.5.4-linux-x86_64.tar.gz gringo-4.5.4-linux-x86_64/gringo --strip-components 1