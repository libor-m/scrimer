Install non python dependencies
===============================
Sequence of commands used to install the non-python dependencies.

Software direcotry
------------------
mkdir sw
cd sw

bedtools
--------
cd ~/sw
wget -O - https://github.com/arq5x/bedtools2/releases/download/v2.22.0/bedtools-2.22.0.tar.gz|tar xvz
cd bedtools2/
make
sudo cp -R bin /opt/bedtools
echo /opt/bedtools/bin >> paths

samtools and tabix
------------------
# samtools and htslib do not support optional installation prefix
# so they're installed in 
cd ~/sw
wget -O - 'http://downloads.sourceforge.net/project/samtools/samtools/1.1/samtools-1.1.tar.bz2?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fsamtools%2Ffiles%2Fsamtools%2F1.1%2F&ts=1416174418&use_mirror=netcologne'|tar xvj
cd samtools-1.1
sudo apt-get install -y ncurses-dev
make && sudo make install

# tabix, bgzip
cd htslib-1.1
make && sudo make install

# vcfutils.pl
wget -O - http://sourceforge.net/projects/samtools/files/samtools/1.1/bcftools-1.1.tar.bz2/download|tar xvj
cd bcftools-1.1/ 
make && sudo make install

lastz
-----
cd ~/sw
echo /opt/lastz/bin >> paths
wget -O - http://www.bx.psu.edu/miller_lab/dist/lastz-1.02.00.tar.gz|tar xvz
cd lastz-distrib-1.02.00/
# remove -Werror from the makefile
nano src/Makefile
make
LASTZ_INSTALL=/opt/lastz/bin sudo -E make install


gmap
----
cd ~/sw
echo /opt/gmap/bin >> paths
wget -O - http://research-pub.gene.com/gmap/src/gmap-gsnap-2014-10-22.tar.gz|tar xvz
cd gmap-2014-10-22/
./configure --prefix=/opt/gmap
make
sudo make install

smalt
-----
cd ~/sw
echo /opt/smalt/bin >> paths
git clone http://git.code.sf.net/p/smalt/code smalt
cd smalt/
./configure --prefix=/opt/smalt
make
sudo make install

parallel
--------
# do not install parallel into /opt, but systemwide
cd ~/sw
wget -O - http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2|tar xj
cd parallel-20141022/
./configure
make
sudo make install

blat and ispcr
--------------
cd ~/sw
echo /opt/kent/bin >> paths
sudo apt-get install unzip libpng-dev
wget http://users.soe.ucsc.edu/~kent/src/blatSrc35.zip
wget http://users.soe.ucsc.edu/~kent/src/isPcr.zip
unzip '*.zip'
mkdir kent_bin
export MACHTYPE
export BINDIR=~/sw/kent_bin
cd blatSrc/
make
cd ../isPcrSrc/
cd isPcr/isPcr/

# remove BINDIR=.. and -Werror
nano ../../inc/common.mk

make
cd ~/sw/kent_bin
sudo mkdir -p /opt/kent/bin
sudo cp * /opt/kent/bin

Primer3
-------
cd ~/sw
echo /opt/primer3/bin >> paths
wget -O - 'http://downloads.sourceforge.net/project/primer3/primer3/2.3.6/primer3-src-2.3.6.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fprimer3%2Ffiles%2F&ts=1416259162&use_mirror=heanet'|tar xvz
cd primer3-2.3.6/src
make
sudo mkdir -p /opt/primer3/bin
sudo cp long_seq_tm_test ntdpal ntthal oligotm primer3_core /opt/primer3/bin
sudo cp -R primer3_config /opt/primer3

cutadapt
--------
pip install cutadapt

FastQC
------
cd ~/sw
echo /opt/FastQC >> paths
sudo apt-get install openjdk-7-jre
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip
unzip fastqc_v0.11.2.zip
sudo cp -R FastQC /opt
sudo chmod +x /opt/FastQC/fastqc
 
agrep and tre-agrep
-------------------
cd ~/sw
git clone https://github.com/Wikinaut/agrep.git
cd agrep/
make
sudo cp agrep /usr/local/bin

sudo apt-get install -y tre-agrep 

sort-alt
--------
cd ~/sw
git clone https://github.com/lh3/foreign.git
cd foreign/sort/
make
mv sort sort-alt
sudo cp sort-alt /usr/local/bin

pv
--
cd ~/sw
wget -O - http://www.ivarch.com/programs/sources/pv-1.5.7.tar.bz2|tar xvj
cd pv-1.5.7/
./configure
make && sudo make install

mawk
----
wget -O - http://invisible-island.net/datafiles/release/mawk.tar.gz|tar xvz
cd mawk-1.3.4-20141027/
./configure
make && sudo make install

bit.ly data hacks
-----------------
. ~/scrimer-env/bin/activate
pip install data_hacks

freebayes
---------
sudo apt-get install cmake
git clone --recursive git://github.com/ekg/freebayes.git
cd freebayes/
make && sudo make install
