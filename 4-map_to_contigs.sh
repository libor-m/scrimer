#------------------------------------------------
# Operation
#------------------------------------------------

#  1. map contigs to scaffold
#  2. merge mapping output to a single file

# tools
TOOLS=/opt/samtools-0.1.18:~/data/sw_testbed/IGVTools:~/data/sw_testbed/smalt-0.7.0.1
export PATH=$TOOLS:$PATH


# 1. map contigs to scaffold
#------------------------------------------------ 

#------------------------------------------------ smalt
SCAFFOLD=33-scaffold/lx4.fasta 
INFILES=12-cutadapt/*.fastq
OUT=40-map-smalt
SMALT_IDX=${SCAFFOLD%/*}/smalt/k13s4

# create smalt index
mkdir -p ${SMALT_IDX%/*}
smalt index -s 4 $SMALT_IDX $SCAFFOLD

# map each file, smalt is multithreaded so feed the files sequentially
mkdir -p $OUT
for FQFILE in $INFILES
do
  SAMFILE=$OUT/$( basename ${FQFILE%.*} ).sam
  smalt map -n 8 -p -f sam -o $SAMFILE $SMALT_IDX $FQFILE
done

#------------------------------------------------ gsnap?

# 2. merge mapping output to a single file 
#------------------------------------------------ 

# create a fasta index
samtools faidx $SCAFFOLD

# according to real sample info, create a readgroups file
# samtools merge -r switch attaches a read group to each alignment (line) in input 
# according to original filename
# so the format is:
# @RG	ID:$BASENAME	SM:$SAMPLE	LB:${BASENAME%%.*}	PL:LS454
# the library name (LB) is important because of rmdup
# description (DS) is used to identify the species
# FIXME: this does not work when copy-pasting to bash, 
# because shell interprets the tab character (autocomplete)
# 
# the order matters for the vcf output,
# the sample columns order is probably the order of first apperance in the @RG
#

# template generator
# reorder lines and fill the '??' marked places

# this won't work with filenames conatining spaces
echo $INFILES|tr " " "\n"|xargs -n1 basename|sed s/.fastq//|gawk '{OFS="\t";print "@RG", "ID:" $0, "SM:??", "LB:" gensub(/\..*$/,"",$0), "PL:LS454", "DS:??";}' > $OUT/readgroups.txt

# alternative using find
DIR=12-cutadapt
find $DIR -name '*.fastq'|xargs -n1 basename|sed s/.fastq//|gawk '{OFS="\t";print "@RG", "ID:" $0, "SM:??", "LB:" gensub(/\..*$/,"",$0), "PL:LS454", "DS:??";}' > $OUT/readgroups.txt

# generate the sequence headers, using any (first) .sam file
export SAMFILE=$( echo $OUT/*.sam|gawk '{print $1;}' )
samtools view -S -t $SCAFFOLD.fai -H $SAMFILE > $OUT/sequences.txt
cat $OUT/sequences.txt $OUT/readgroups.txt > $OUT/sam-header.txt

# samtools merge requires sorted alignments, sort them in parallel
parallel "samtools view -but $SCAFFOLD.fai {} | samtools sort - {.}" ::: $OUT/*.sam

# merge all the alignments, not removing the dupes?
# http://seqanswers.com/forums/showthread.php?t=6543
# http://seqanswers.com/forums/showthread.php?t=5424
# using /[GH]*.bam to avoid generated files (like alldup.bam) in the expansion
samtools merge -ru -h $OUT/sam-header.txt - $OUT/[GH]*.bam | samtools sort - $OUT/alldup
samtools index $OUT/alldup.bam

#------------------------------------------------
# Visualization
#------------------------------------------------

# unmapped read counts
parallel 'echo $( cut -f2 {}|grep -c "^4$" ) {}' ::: $OUT/*.sam

# mapping statistics
samtools idxstats $OUT/alldup.bam|gawk '{map += $3; unmap += $4;} END {print  unmap/map;}'

# coverage sums for IGV
igvtools count -z 5 -w 25 -e 250 $OUT/alldup.bam  $OUT/alldup.bam.tdf ${CONTIGS%.*}.genome

#------------------------------------------------
# spare parts
#------------------------------------------------

# header fix
samtools reheader $OUT/sam-header2.txt $OUT/alldup.bam > $OUT/alldup2.bam

# count the coverage for wider zoom levels in igv 
igvtools count -z 5 -w 25 -e 250 "${SAMFILE%.*}".bam  "${SAMFILE%.*}".bam.tdf ${CONTIGS%.*}.genome

##### single file process 
# easiest way to get a correct sequences header -> choose the smallest .sam
# (also useful for testing options)

# convert, sort and index the sam
SAMFILE=$OUT/$( basename ${FQFILE%.*} ).sam
samtools view -but $SCAFFOLD.fai $SAMFILE | samtools sort - "${SAMFILE%.*}"
samtools index "${SAMFILE%.*}".bam

# goby test, takes 1/5 of .bam (220 MB to 1.1 GB), with --preserve* there is no space gain
# IGV cannot use the sample info form goby to separate the groups regretably..
samtools fillmd -u $OUT/alldup.bam $SCAFFOLD| goby 4g stc -i - -o $OUT/alldup.goby -g $SCAFFOLD --preserve-all-tags --preserve-all-mapped-qualities --sorted

# sort goby if needed
goby 4g sort -t -1 -o $OUT/alldup.goby.sort $OUT/alldup.goby
