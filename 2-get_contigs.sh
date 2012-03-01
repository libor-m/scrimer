# yasra requires single fasta file as input
# fastx toolkit try -- didn't like quality of '0'
cat 14-cutadapt-2/*.fastq | fastq_to_fasta -n > 24-yasra-zfmrna/allreads.fa

# fastx failed, by-hand-solution
cat 14-cutadapt-2/*.fastq | awk '(NR % 4) == 1 {sub(/^@/, ">", $1); print $1;} (NR % 4) == 2 { print; }' > 24-yasra-zfmrna/allreads.fa

# run yasra -- the Makefile had to be edited to include [unmask] action, 
# because mrna/est reference files are in lowercase
cd 24-yasra-zfmrna && make single_step > Summary

samtools faidx 24-yasra-zfmrna/Final_Assembly
SAMFILE=24-yasra-zfmrna/alignments.sam
samtools view -but 24-yasra-zfmrna/Final_Assembly.fai $SAMFILE|samtools sort - ${SAMFILE%%.*}.bam
samtools index 24-yasra-zfmrna/alignments.bam

# IGV gui must be used to create a '.genome' via 'Import Genome ...'
# http://www.broadinstitute.org/igv/LoadGenome
~/data/sw_testbed/IGVTools/igvtools count -z 5 -w 25 -e 250 "${SAMFILE%%.*}".bam  "${SAMFILE%%.*}".bam.tdf taeGut1

# check the self similarity of the contigs
# using lastz
SEQFILE=0a-jp-newbler-contigs/lu_master500.fasta
# using [multiple] breaks the --self, --nomirror, --rdotplot options
lastz $SEQFILE[multiple] $SEQFILE --nomirror --gfextend --chain --gapped --entropy --format=sam --rdotplot=${SEQFILE%%.*}_self.tsv --progress> ${SEQFILE%%.*}_self.sam

# better like this, --entropy to downscore polyA/T runs, softsam to keep the sequence context for general overview of the hit
# or use CIGAR?
lastz $SEQFILE[multiple] $SEQFILE --gfextend --chain --gapped --entropy --format=softsam --progress=1000> ${SEQFILE%%.*}_self.sam

samtools view -buS $SAMFILE|samtools sort - "${SAMFILE%%.*}"
samtools index "${SAMFILE%%.*}.bam"

# the lastz command does not remove the 'trivial' matches, because --self does not work with [multiple]
# it's necessary to remove those by hand
gawk -F "\t" '(NF < 10 || (NF > 10 && $1 != $3))' $SAMFILE|samtools view -buS -|samtools sort - "${SAMFILE%%.*}_f"
samtools index "${SAMFILE%%.*}_f.bam"

