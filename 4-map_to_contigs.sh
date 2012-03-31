# data
CONTIGS=0a-jp-newbler-contigs/lu_master500.fasta
OUTDIR=40-map-smalt

# tools
TOOLS=/opt/samtools-0.1.18:~/data/sw_testbed/IGVTools:
export PATH=$TOOLS:$PATH

######
# use smalt to map the reads to the contigs
######

# build the index in the directory of the genome
# recommended settings for 454 (step 4, k-mer size 13)
mkdir -p ${CONTIGS%/*}/smalt
smalt index -s 4 ${CONTIGS%/*}/smalt/k13s4 $CONTIGS

# smalt seems to be a lot less effective when there is a lot of sequences in the 'genome'
# probably redo the alignments with some coverage threshold, to avoid 200H 20M 300H matches
# map the reads in a loop
mkdir -p $OUTDIR
for FQFILE in 14-cutadapt-2/*.fastq
do
  SAMFILE=$OUTDIR/$( basename ${FQFILE%.*} ).sam
  smalt map -n 8 -p -f sam -o $SAMFILE ${CONTIGS%/*}/smalt/k13s4 $FQFILE
done

# all/mapped read counts
parallel 'echo $( cut -f2 {}|grep -c "^4$" ) {}' ::: $OUTDIR/*.sam


######
# samtools chiseling of the output
######

# create a fasta index
samtools faidx $CONTIGS

##### single file process 
# easiest way to get a correct sequences header -> choose the smallest .sam
# (also useful for testing options)

# convert, sort and index the sam
SAMFILE=$OUTDIR/$( basename ${FQFILE%.*} ).sam
samtools view -but $CONTIGS.fai $SAMFILE | samtools sort - "${SAMFILE%.*}"
samtools index "${SAMFILE%.*}".bam

# count the coverage for wider zoom levels in igv 
# (pseudogenome would be nice for this, because the -e 250 overshoots each read and reports errors
# and the tdf is soo huge..)
igvtools count -z 5 -w 25 -e 250 "${SAMFILE%.*}".bam  "${SAMFILE%.*}".bam.tdf ${CONTIGS%.*}.genome

# generate the sequence headers, using any .sam file
samtools view -S -t $CONTIGS.fai -H $SAMFILE > $OUTDIR/sequences.txt

# according to real sample info, create a readgroups file
# samtools merge -r switch attaches a read group to each alignment (line) in input 
# according to original filename
# so the format is:
# @RG	ID:$BASENAME	SM:$SAMPLE	LB:${BASENAME%%.*}	PL:LS454
# the library name (LB) is important because of rmdup
# description (DS) is used to identify the species
# FIXME: this does not work when copy-pasting, 
# because shell interprets the tab character (autocomplete)
# 
# the order matters for the vcf output,
# the sample columns order is probably the order of first apperance in the @RG
#
cat > $OUTDIR/readgroups.txt <<'EOF'
@RG	ID:G59B7NP01	SM:lu02	LB:G59B7NP01	PL:LS454	DS:ll
@RG	ID:G60Z2EH01.RL9	SM:lu04	LB:G60Z2EH01	PL:LS454	DS:lm
@RG	ID:G60Z2EH01.RL11	SM:lu06	LB:G60Z2EH01	PL:LS454	DS:lm
@RG	ID:G60Z2EH02.RL10	SM:lu05	LB:G60Z2EH02	PL:LS454	DS:ll
@RG	ID:G60Z2EH02.RL12	SM:lu07	LB:G60Z2EH02	PL:LS454	DS:ll
@RG	ID:GS60IET02.RL1	SM:lu01	LB:GS60IET02	PL:LS454	DS:lm
@RG	ID:GS60IET02.RL2	SM:lu03	LB:GS60IET02	PL:LS454	DS:lm
@RG	ID:GSVZDOM01.RL1	SM:lu01	LB:GSVZDOM01	PL:LS454	DS:lm
@RG	ID:GSVZDOM01.RL2	SM:lu03	LB:GSVZDOM01	PL:LS454	DS:lm
@RG	ID:GSVZDOM02	SM:lu02	LB:GSVZDOM02	PL:LS454	DS:ll
@RG	ID:GYAB93P02.RL1	SM:lu01	LB:GYAB93P02	PL:LS454	DS:lm
@RG	ID:GYAB93P02.RL2	SM:lu03	LB:GYAB93P02	PL:LS454	DS:lm
EOF

# samtools merge requires sorted alignments, sort them in parallel
parallel "samtools view -but $CONTIGS.fai {} | samtools sort - {.}" ::: $OUTDIR/*.sam

# generate the sequence headers, using any .sam file
samtools view -S -t $CONTIGS.fai -H $SAMFILE > $OUTDIR/sequences.txt
cat $OUTDIR/sequences.txt readgroups.txt > $OUTDIR/sam-header.txt

# merge all the alignments, not removing the dupes?
# http://seqanswers.com/forums/showthread.php?t=6543
# http://seqanswers.com/forums/showthread.php?t=5424
# using /G*.bam to avoid generated files (like alldup.bam) in the expansion
samtools merge -ru -h $OUTDIR/sam-header.txt - $OUTDIR/G*.bam | samtools sort - $OUTDIR/alldup
samtools index $OUTDIR/alldup.bam

# final mapping statistics
samtools idxstats $OUTDIR/alldup.bam|gawk '{map += $3; unmap += $4;} END {print  unmap/map;}'

# coverage sums for IGV
igvtools count -z 5 -w 25 -e 250 $OUTDIR/alldup.bam  $OUTDIR/alldup.bam.tdf ${CONTIGS%.*}.genome

# header fix
samtools reheader $OUTDIR/sam-header2.txt $OUTDIR/alldup.bam > $OUTDIR/alldup2.bam
