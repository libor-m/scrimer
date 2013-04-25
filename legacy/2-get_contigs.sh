#------------------------------------------------
# Operation
#------------------------------------------------

#  1. use newbler to assemble the files
#  2. remove contigs that are similar to each other

# tools used: newbler, lastz

#  1. use newbler to assemble the files
#------------------------------------------------

# data from previous steps
SEQFILE=$OUT/lu_master500_v2.fna

OUT=20-jp-contigs

#TODO: add some working newbler command 
# - current assembly was done by jpaces

# get one exclusive zewura, interactive from screen at skirit
qsub -q q_2d@wagap.cerit-sc.cz -l mem=400gb -l nodes=1:ppn=80:nodecpus80:cl_zewura#excl -I

# add newbler to path
export PATH=~/brno3/newbler-2.6/bin:$PATH

# run test assembly 
runAssembly -o 21-newbler-test -cdna -cpu 19 12-cutadapt/GSVZDOM01.RL1.fastq

# run full assembly 
runAssembly -o 22-newbler -cdna -cpu 19 -m 12-cutadapt/*.fastq

#  2. remove contigs that are similar to each other
#------------------------------------------------

# taken from lastz human-chimp example, should be report only highly similar hits
# filter out self matches with awk
# takes 8 minutes, finds 190K pairs
lastz $SEQFILE[multiple] $SEQFILE \
      --step=10 --seed=match12 --notransition --exact=20 --noytrim \
      --match=1,5 --ambiguous=n \
      --coverage=90 --identity=95 \
      --format=general:name1,size1,start1,name2,size2,start2,strand2,identity,coverage \
      | awk '($1 != $4)' > $SEQFILE.lastz-self


# find the redundant sequences
tail -n +2 $SEQFILE.lastz-self|./find_redundant_sequences.py > $SEQFILE.redundant

# add the short ones
grep '>' $SEQFILE|awk '{ sub(/length=/,"",$3); sub(/^>/, "", $1); if($3 < 300) print $1;}' >> $SEQFILE.redundant

# get rid of the redundant ones
./seq_filter_by_id.py $SEQFILE.redundant 1 $SEQFILE fasta - $SEQFILE.filtered

#------------------------------------------------
# Visualization
#------------------------------------------------

# check the contig length distribution
grep '>' $SEQFILE|awk '{ sub(/length=/,"",$3); print $3}' > $SEQFILE.lengths

# view how many contigs we got after filtering
grep -c '>' $SEQFILE.filtered

#------------------------------------------------
# spare parts
#------------------------------------------------

# using [multiple] breaks the --self, --nomirror, --rdotplot options
lastz $SEQFILE[multiple] $SEQFILE --nomirror --gfextend --chain --gapped --entropy --format=sam --rdotplot=${SEQFILE%%.*}_self.tsv --progress> ${SEQFILE%%.*}_self.sam

# better like this, --entropy to downscore polyA/T runs, softsam to keep the sequence context for general overview of the hit
# or use CIGAR?
lastz $SEQFILE[multiple] $SEQFILE --gfextend --chain --gapped --entropy --format=softsam --progress=1000> ${SEQFILE%%.*}_self.sam


# IGV gui must be used to create a '.genome' via 'Import Genome ...'
# http://www.broadinstitute.org/igv/LoadGenome
~/data/sw_testbed/IGVTools/igvtools count -z 5 -w 25 -e 250 "${SAMFILE%%.*}".bam  "${SAMFILE%%.*}".bam.tdf taeGut1


samtools view -buS $SAMFILE|samtools sort - "${SAMFILE%%.*}"
samtools index "${SAMFILE%%.*}.bam"

# the lastz command does not remove the 'trivial' matches, because --self does not work with [multiple]
# it's necessary to remove those by hand
gawk -F "\t" '(NF < 10 || (NF > 10 && $1 != $3))' $SAMFILE|samtools view -buS -|samtools sort - "${SAMFILE%%.*}_f"
samtools index "${SAMFILE%%.*}_f.bam"
