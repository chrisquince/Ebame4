#!/bin/bash

#Install megahit
cd ~
mkdir bin
mkdir repos
cd ~/repos
git clone https://github.com/voutcn/megahit
cd megahit/
sudo apt-get install --assume-yes zlib1g-dev
make
cp megahit* ~/bin

# Install BWA
cd ~/repos
git clone https://github.com/lh3/bwa.git
cd bwa; make
cp bwa ~/bin

#bam-readcount
cd ~/repos
sudo apt-get install --assume-yes build-essential git-core cmake zlib1g-dev libncurses-dev patch
git clone https://github.com/genome/bam-readcount.git
mkdir bam-readcount-build
cd bam-readcount-build/
cmake ../bam-readcount
make
cp bin/bam-readcount ~/bin/

#samtools
cd ~/repos
wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
tar xvfj samtools-1.3.1.tar.bz2 
cd samtools-1.3.1/ 
sudo apt-get install --assume-yes libcurl4-openssl-dev libssl-dev
./configure --enable-plugins --enable-libcurl --with-plugin-path=$PWD/htslib-1.3.1
make all plugins-htslib
cp samtools ~/bin/  

#bedtools
cd ~/repos
sudo apt-get install --assume-yes bedtools

#prodigal
cd ~/repos
wget https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux 
cp prodigal.linux ~/bin/prodigal
chmod +rwx ~/bin/prodigal

#gnu parallel
sudo apt-get install --assume-yes parallel
sudo apt install --assume-yes fasttree
 
#mafft
sudo apt install --assume-yes mafft

#standalone blast
cd ~/repos
wget https://desmandatabases.s3.climb.ac.uk/ncbi-blast-2.5.0+-x64-linux.tar.gz

tar -xvzf ncbi-blast-2.5.0+-x64-linux.tar.gz

cp ncbi-blast-2.5.0+/bin/* ~/bin

#diamond
cd ~/repos
mkdir diamond
cd diamond
wget http://github.com/bbuchfink/diamond/releases/download/v0.9.22/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
cp diamond ~/bin

#install R
cd ~/repos
sudo echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | sudo tee -a /etc/apt/sources.list
gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
 gpg -a --export E084DAB9 | sudo apt-key add -
 sudo apt-get update
sudo apt-get install --assume-yes r-base r-base-dev


cd ~/repos/
git clone https://github.com/scapella/trimal
cd trimal/source
make
cp trimal ~/bin/

#install Spades
cd ~
mkdir Storage
sudo mount /dev/vdb Storage
sudo chown -R ubuntu:ubuntu Storage

mkdir ~/Storage/Installation
cd ~/Storage/Installation
wget http://cab.spbu.ru/files/release3.10.1/SPAdes-3.10.1-Linux.tar.gz
tar -xzf SPAdes-3.10.1-Linux.tar.gz
cd SPAdes-3.10.1-Linux/bin/
echo 'export PATH=~/Storage/Installation/SPAdes-3.10.1-Linux/bin:$PATH' >> ~/.bashrc 

#Install databases
mkdir ~/Storage/Databases

cd ~/Storage/Databases
wget https://desmandatabases.s3.climb.ac.uk/rpsblast_cog_db.tar.gz
tar -xvzf rpsblast_cog_db.tar.gz


#Install NR
#cd ~/Storage/Databases
#mkdir NR
#cd NR
#wget https://desmandatabases.s3.climb.ac.uk/nr.faa
#diamond makedb --in nr.faa -d nr

#wget https://desmandatabases.s3.climb.ac.uk/gi_taxid_prot.dmp
#wget https://desmandatabases.s3.climb.ac.uk/all_taxa_lineage_notnone.tsv

cd ~/Storage/Databases
wget https://kegg.s3.climb.ac.uk/KeggUpdate.tar.gz
tar -xvzf KeggUpdate.tar.gz

#Stopped here

cd ~/repos
sudo apt-get -y install python-pip python3-pip 
pip install cython numpy scipy biopython pandas pip scikit-learn pysam bcbio-gff
pip3 install cython numpy scipy biopython pandas pip scikit-learn pysam bcbio-gff

#install gsl
sudo apt-get install --assume-yes libgsl-dev

#install CONCOCT development branch
cd ~/repos
git clone https://github.com/BinPro/CONCOCT.git
cd CONCOCT
git fetch
git checkout develop
sudo python ./setup.py install

#Install DESMAN
cd ~/repos
git clone https://github.com/chrisquince/DESMAN.git
cd DESMAN
sudo python3 ./setup.py install
cd ..
git clone https://github.com/chrisquince/Ebame4.git
git clone https://github.com/chrisquince/MAGAnalysis.git

#Added these lines to .bashrc
echo 'export CONCOCT=~/repos/CONCOCT' >> ~/.bashrc 
echo 'export DESMAN=~/repos/DESMAN' >> ~/.bashrc 
echo 'export PATH=~/repos/Ebame4/scripts:$PATH' >> ~/.bashrc 


cd ~/Storage/
  
mkdir Data
cd Data
mkdir InfantGut
cd InfantGut
wget https://infantgut.s3.climb.ac.uk/ReadsSub.tar.gz
tar -xvzf ReadsSub.tar.gz 

cd ~/Storage/
mkdir Projects

cd ~
ln -s Storage/Projects .

cd ~/Projects


rm ~/Storage/Data/InfantGut/ReadsSub/sample8_R*.fastq
rm ~/Storage/Data/InfantGut/ReadsSub/sample9_R*.fastq
rm ~/Storage/Data/InfantGut/ReadsSub/sample10_R*.fastq
rm ~/Storage/Data/InfantGut/ReadsSub/sample11_R*.fastq

mkdir InfantGut
ln -s ~/Storage/Data/InfantGut/ReadsSub ~/Projects/InfantGut

#Download minikraken database
cd ~/Storage/Databases/
wget https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_8GB.tgz
tar -xzvf minikraken_20171019_8GB.tgz 

#Install Kraken
cd ~/repos/
git clone https://github.com/DerrickWood/kraken.git
cd kraken
./install_kraken.sh ~/bin/

cd ~

sudo apt install --assume-yes mercurial

pip install matplotlib --user

cd ~/Installation

curl -O https://raw.githubusercontent.com/IFB-ElixirFr/biosphere-commons/master/scripts/utils/biosphere-nextcloud-client.sh

source ./biosphere-nextcloud-client.sh

sudo ~/repos/Ebame4/scripts/RPInstall.sh ggplot2
sudo ~/repos/Ebame4/scripts/RPInstall.sh reshape
sudo ~/repos/Ebame4/scripts/RPInstall.sh reshape2
sudo ~/repos/Ebame4/scripts/RPInstall.sh gplots
sudo ~/repos/Ebame4/scripts/RPInstall.sh getopt
sudo ~/repos/Ebame4/scripts/RPInstall.sh vegan
sudo ~/repos/Ebame4/scripts/RPInstall.sh ellipse
sudo ~/repos/Ebame4/scripts/RPInstall.sh plyr
sudo ~/repos/Ebame4/scripts/RPInstall.sh grid
sudo ~/repos/Ebame4/scripts/RPInstall.sh gridExtra

sudo apt install --assume-yes mafft
cd ~/Storage/Databases
wget https://kegg.s3.climb.ac.uk/Combined.tar.gz 
tar -xvzf Combined.tar.gz 
