#
# obsolote code, just for reference, use the #-step.sh scripts
#
######
# use ssahaEST to perform the mapping to zebra finch genome
######
# necessary for the program to locate MTX dir, can be used to 
# issue shorter commands  
export SSAHADIR_EST=~/data/sw_testbed/ssahaEST_1_0_1_b.amd64
# error in the readme file..
export SSAHAEST_DIR=~/data/sw_testbed/ssahaEST_1_0_1_b.amd64

# first issued command generates the hash index
$SSAHAEST_DIR/ssahaEST $FQFILE /data/genomes/taeGut1/taeGut1.fa -align 0 -output cigar -save /data/genomes/taeGut1/ssahaEST/default > ...

export CIGFILE=20-map-ssahaEST/$( basename ${FQFILE%%.*} ).cigar
$SSAHAEST_DIR/ssahaEST $FQFILE /data/genomes/taeGut1/taeGut1.fa -align 0 -output cigar > $CIGFILE

# map the reads, really... this segfaults on xukol
$SSAHAEST_DIR/ssahaEST $FQFILE /data/genomes/taeGut1/taeGut1.fa -align 0 -output cigar -use /data/genomes/taeGut1/ssahaEST/default 

# sshaEST is just broken and does not seem to have much support on the web

######
# use smalt to map the reads
######
export SMALT_DIR=~/data/sw_testbed/smalt-0.5.8

# build the index, recommended settings for 454 (step 4, k-mer size 13)
$SMALT_DIR/smalt index -s 4 /data/genomes/taeGut1/smalt/taeGut1k13s4 /data/genomes/taeGut1/taeGut1.fa

# map the reads
export SAMFILE=21-map-smalt/$( basename ${FQFILE%.*} ).sam
$SMALT_DIR/smalt map -n 8 -f sam -o $SAMFILE /data/genomes/taeGut1/smalt/taeGut1k13s4 $FQFILE

# map the reads in a loop
for FQFILE in 13-cutadapt/*.fastq
do
  SAMFILE=23-map-smalt-p/$( basename ${FQFILE%.*} ).sam
  $SMALT_DIR/smalt map -n 8 -p -f sam -o $SAMFILE /data/genomes/taeGut1/smalt/taeGut1k13s4 $FQFILE
done

# all/mapped reads counts
wc -l $SAMFILE && cut -f2 $SAMFILE|grep -c "^4$"

parallel 'echo $( cut -f2 {}|grep -c "^4$" ) {}' ::: 23-map-smalt-p/*.sam

######
# samtools chiseling of the output
######
# to view the alignments in igv, BAM is recommended and must be indexed by igvtools
export PATH=/opt/samtools-0.1.18:$PATH

# create a fasta index
samtools faidx /data/genomes/taeGut1/taeGut1.fa

# convert, sort and index the sam
samtools view -but /data/genomes/taeGut1/taeGut1.fa.fai $SAMFILE | samtools rmdup -s - -| samtools sort - "${SAMFILE%%.*}"
samtools index "${SAMFILE%%.*}".bam

# count the coverage for high zoom levels in igv
~/data/sw_testbed/IGVTools/igvtools count -z 5 -w 25 -e 250 "${SAMFILE%%.*}".bam  "${SAMFILE%%.*}".bam.tdf taeGut1

######

# parallel processing 
parallel "samtools view -but /data/genomes/taeGut1/taeGut1.fa.fai {} | samtools sort - {.}" ::: 23-map-smalt-p/*.sam

# generate the headers, using a dummy .sam file
samtools view -S -t /data/genomes/taeGut1/taeGut1.fa.fai -H 23-map-smalt-p/GSVZDOM01.RL1.sam > 23-map-smalt-p/sq.txt
cat  23-map-smalt-p/sq.txt 23-map-smalt-p/rg_ll.txt > 23-map-smalt-p/ll_hdr.txt
cat  23-map-smalt-p/sq.txt 23-map-smalt-p/rg_lm.txt > 23-map-smalt-p/lm_hdr.txt
cat  23-map-smalt-p/sq.txt 23-map-smalt-p/rg_all.txt > 23-map-smalt-p/all_hdr.txt

# merge per species
cd 23-map-smalt-p
samtools merge -ru -h ll_hdr.txt - G59B7NP01.bam G60Z2EH02.RL10.bam G60Z2EH02.RL12.bam GSVZDOM02.bam | samtools rmdup -s - -| samtools sort - ll
# [bam_rmdupse_core] 580280 / 1427457 = 0.4065 in library 'll'
samtools index ll.bam

# whole python script adding .bam extension is still shorter than the list;)
LMIDS="G60Z2EH01.RL9 G60Z2EH01.RL11 GS60IET02.RL1 GS60IET02.RL2 GSVZDOM01.RL1 GSVZDOM01.RL2 GYAB93P02.RL1 GYAB93P02.RL2"
samtools merge -ru -h lm_hdr.txt - $( python -c "print ' '.join([x + '.bam' for x in '$LMIDS'.split()])" ) | samtools rmdup -s - -| samtools sort - lm
# [bam_rmdupse_core] 737148 / 1782684 = 0.4135 in library 'lm'
samtools merge -ru -h lm_hdr.txt - $( python -c "print ' '.join([x + '.bam' for x in '$LMIDS'.split()])" ) | samtools sort - lmdup

# all in one
samtools merge -ru -h all_hdr.txt - G*.bam | samtools rmdup -s - -| samtools sort - all
# samtools rmdup works according to @RG LB: -- take care
samtools merge -ru -h all_hdr.txt - G*.bam | samtools sort - alldup

# final mapping statistics
samtools idxstats 23-map-smalt-p/alldup.bam|gawk '{map += $3; unmap += $4;} END {print  unmap/map;}'