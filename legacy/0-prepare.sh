#------------------------------------------------
# Operation
#------------------------------------------------

#  1. downlaod and prepare the reference genome
#  2. downlaod all needed annotations
#  3. and build genome indexes for all programs used in the pipeline
#  4. quality check of the input data

# tools used: bedtools, tabix, sort-alt, gmap, samtools, smalt, fastqc

# 1. reference genome
#------------------------------------------------

# database list at http://hgdownload.cse.ucsc.edu/downloads.html#zebrafinch
# download full data set, interrupt during xenoMrna (not needed, too big)
# command displayed in the file listing

# location of genome data that can be sharad among users
cd /data/genomes
GENOME=taeGut1
mkdir $GENOME
cd $GENOME
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/$GENOME/bigZips/ .

# check downloaded data integrity
md5sum -c md5sum.txt
cat *.md5|md5sum -c

# unpack the genome - this differs for genomes
# some are in single .fa, some are split by chromosomes
tar xvzf chromFa.tar.gz

# create concatenated chicken genome, use Heng Li's sort-alt
# to get common ordering of chromosomes
find chromFa -type f|sort-alt -N|xargs cat > galGal3.fa

# 2. annotataions
#------------------------------------------------
# annotation data is best obtained in UCSC table browser
# in BED format and then sorted and indexed by BEDtools
# http://genome.ucsc.edu/cgi-bin/hgTables?db=taeGut1
ANNOT=annot
sortBed -i $ANNOT/ensGene.bed > $ANNOT/ensGene.sorted.bed
bgzip $ANNOT/ensGene.sorted.bed
tabix -p bed $ANNOT/ensGene.sorted.bed.gz

# 
# 3. indexes
#------------------------------------------------
#TODO: more commands

# index chromosome positions in the genome file for samtools
samtools faidx /data/genomes/taeGut1/taeGut1.fa

# build gmap index for zebra finch
gmap_build -d gmap_taeGut1 -D /data/genomes/taeGut1 /data/genomes/taeGut1/taeGut1.fa


# smalt index
# recommended settings for 454 (step 4, k-mer size 13)
/data/genomes/taeGut1/smalt/taeGut1k13s4
mkdir -p smalt
smalt index -s 4 smalt/${GENOME}k13s4 $GENOME.fa

# 4. quality check
#------------------------------------------------
mkdir 01-fastqc
fastqc --outdir=01-fastqc --noextract --threads=8 00-raw/*.fastq
