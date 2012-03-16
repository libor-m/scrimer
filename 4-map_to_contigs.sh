# dirs
CONTIGS=0a-jp-newbler-contigs/lu_master500.fasta
OUTDIR=40-map-smalt

# tools
TOOLS=/opt/samtools-0.1.18:~/data/sw_testbed/IGVTools
export PATH=$TOOLS:$PATH

######
# use smalt to map the reads to the contigs
######

# build the index, recommended settings for 454 (step 4, k-mer size 13)
smalt index -s 4 $OUTDIR/lu_master500/k13s4 $CONTIGS

# smalt seems to be a lot less effective when there is a lot of sequences in the 'genome'
# map the reads in a loop
for FQFILE in 14-cutadapt-2/*.fastq
do
  SAMFILE=$OUTDIR/$( basename ${FQFILE%.*} ).sam
  smalt map -n 8 -p -f sam -o $SAMFILE 40-map-smalt/lu_master500/k13s4 $FQFILE
done

# all/mapped read counts
parallel 'echo $( cut -f2 {}|grep -c "^4$" ) {}' ::: $OUTDIR/*.sam


######
# samtools chiseling of the output
######

# create a fasta index
samtools faidx $CONTIGS

##### single file process (useful for testing options)

# convert, sort and index the sam
SAMFILE=$OUTDIR/$( basename ${FQFILE%.*} ).sam
samtools view -but $CONTIGS.fai $SAMFILE | samtools rmdup -s - -| samtools sort - "${SAMFILE%.*}"
samtools index "${SAMFILE%.*}".bam

# count the coverage for wider zoom levels in igv 
# (pseudogenome would be nice for this, because the -e 250 overshoots each read and reports errors, and the tdf is huge..)
igvtools count -z 5 -w 25 -e 250 "${SAMFILE%.*}".bam  "${SAMFILE%.*}".bam.tdf ${CONTIGS%.*}.genome


###### parallel processing 
parallel "samtools view -but $CONTIGS.fai {} | samtools sort - {.}" ::: 23-map-smalt-p/*.sam

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