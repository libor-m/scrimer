#------------------------------------------------
# Operation
#------------------------------------------------

#  1. map contigs to reference genome (gmap, sim4db, exonerate)
#  2. transfer genome annotations to our contigs
#  3. create 'transcript scaffold' using the annotations

# tools used: samtools, bedtools, tabix, gmap, sim4db, 

# 1. map contigs to reference genome (gmap, sim4db, exonerate)
#------------------------------------------------ 

#------------------------------------------------ gmap
# data from previous steps
INFILE=20-jp-contigs/lu_master500_v2.fna.filtered
GMAP_IDX_DIR=/data/genomes/taeGut1
GMAP_IDX=gmap_taeGut1

OUT=30-tg-gmap
OUTFILE=$OUT/lu_master300_v2.gmap.gff3

# map the jpaces contigs to the zebra finch genome
gmap -D $GMAP_IDX_DIR -d $GMAP_IDX -f gff3_gene -B 3 -x 30 -t 8\
    --cross-species $INFILE  > $OUTFILE

#------------------------------------------------ sim4db
# data from previous steps
INFILE=20-jp-contigs/lu_master500_v2.fna.filtered
GENOME=/data/genomes/taeGut1/taeGut1.fa
SMALT_IDX=/data/genomes/taeGut1/smalt/taeGut1k13s4

OUT=31-tg-sim4db
mkdir $OUT

# these values are derived, it's not necessary to change them
FRAGS=$OUT/${INFILE##*/}.frags
SMALT_OUT=$FRAGS.cigar
SIM4_SCR=${FRAGS%.*}.sim4scr
OUT0=${FRAGS%.*}.tmp.gff3
OUTFILE=${FRAGS%.*}.gff3

# use a fast mapper to find all +-50 KB windows for predicting exon/gene models with sim4db
# smalt with -d on shattered contigs could work reasonably well and fast

# create fragments, using slightly modified fasta_fragments.py from lastz distribution
cat $INFILE | ./fasta_fragments.py --step=80 > $FRAGS

# map the fragments with smalt (takes few minutes), reporting all hits (-d -1) scoring over 60
smalt map -n 8 -f cigar -o $SMALT_OUT -d -1 -m 60 $SMALT_IDX $FRAGS

# construct the script for sim4db
cat $SMALT_OUT|./cigar_to_sim4db_scr.py $GENOME.fai | sort --key=5n,5 > $SIM4_SCR

# run sim4db using the script (takes several seconds for the whole genome !!)
sim4db -genomic $GENOME -cdna $INFILE -script $SIM4_SCR -output $OUT0 -gff3 -interspecies -mincoverage 70 -minidentity 90 -minlength 60 -alignments -threads 7

# fix chromosome names 
sed s/^[0-9][0-9]*:chr/chr/ $OUT0 > $OUTFILE

#------------------------------------------------ exonerate
# TODO: try exonerate, when we get to methods for comparison of generated mappings

# 2. transfer genome annotations to our contigs
#------------------------------------------------ 
# liftover using the exon mapper coordinates
#
# sim4db manual (http://sourceforge.net/apps/mediawiki/kmer/index.php?title=Getting_Started_with_Sim4db)
# Exon coordinates are nucleotide based, starting from 1. Genomic coordinates are always in the original sequence, 
# while the cDNA coordinates will refer to positions in the reverse complement of the sequence if the match orientation is indicated as 'complement'.
#
# --> this is unnecessary, because the orientation of the transcript can be deduced from the target chromosome strand ..?
# --> patched in sim4db, patch on sourceforge
# TODO: add source to github?

# use liftover to transfer the ensGenes exon annotations
# do not trust the mappings per se, they can contain introns 
# - only trust ensGenes or even RefSeq if we need strict conditions
# each contig mapping to genome creates different coordinate system

# add bedtools and samtools to path
TOOLS=/opt/bedtools/bin:/opt/samtools-0.1.18
export PATH=$TOOLS:$PATH

# data from previous steps
# multiple coordinate systems if needed (one system per mapping)
COORDS="30-tg-gmap/lu_master300_v2.gmap.gff3 31-tg-sim4db/lu_master500_v2.fna.filtered.gff3"
# multiple annotations if needed, they're all merged to single gff
ANNOTS=/data/genomes/taeGut1/annot/ensGene_s.bed.gz

# outputs
OUT=32-liftover
mkdir $OUT

for C in $COORDS
do
  ./liftover.py "$C" $ANNOTS > $OUT/${C##*/}-lo.gff3
done  


# 3. create 'transcript scaffold' using the annotations
#------------------------------------------------ 
# construct a 'transcript scaffold' (contigs joined in order of appearance on reference genome chromosomes)
# 'N' gaps should be larger than max read size to avoid the mapping of the reads across gaps

# data from previous steps
INFILE=20-jp-contigs/lu_master500_v2.fna.filtered
ANNOTS=32-liftover/*-lo.gff3

OUT=33-scaffold
mkdir $OUT
GNAME=lx4
OUTGFF=$OUT/$GNAME.gff3

./scaffold.py $INFILE $ANNOTS $OUT/$GNAME.fasta $OUTGFF

# sort, compress and index the merged annotations
# so they can be used further down in the pipeline
OUTFILE=${OUTGFF%.*}.sorted.gff3

sortBed -i $OUTGFF > $OUTFILE
bgzip $OUTFILE
tabix -p gff $OUTFILE.gz


#------------------------------------------------
# Visualization
#------------------------------------------------

# sorted annotations can be opened in IGV

#------------------------------------------------
# spare parts
#------------------------------------------------
# gmap on tophat_test_data
gmap_build -d test_ref -D . -k 12 test_ref.fa
gmap -D . -d test_ref -f samse reads_1.fq reads_2.fq > gmap.sam

# testing various gff formats (BEWARE, IGV depends on the extension being .gff to load the file)
gmap -D . -d test_ref -f gff3_gene test_mrna.fa > gmap.gff3_gene.gff
gmap -D . -d test_ref -f gff3_match_cdna test_mrna.fa > gmap.gff3_match_cdna.gff
gmap -D . -d test_ref -f gff3_match_est test_mrna.fa > gmap.gff3_match_est.gff

# compared to sim4db format
sim4db -genomic test_ref.fa -cdna test_mrna.fa -output sim4db.gff -gff3
# the best option (most similar to sim4db) for gmap is -f gff3_gene

# samtools pipeline for simple aligner testing on tophat_test_data
samtools view -but test_ref.fa.fai gmap.sam|samtools sort -o - tmp|tee gmap.bam|samtools view -h -> gmap_s.sam
samtool index gmap.bam

#------------------------------------------------

gmap -D /data/genomes/taeGut1 -d gmap_taeGut1 -f gff3_gene -B 3 -x 30 -t 6\
    --cross-species $INFILE  > $SAMFILE

#------------------------------------------------
samtools view -but /data/genomes/taeGut1/taeGut1.fa.fai $SAMFILE|samtools sort - "${SAMFILE%%.*}"
samtools index "${SAMFILE%%.*}.bam"
~/data/sw_testbed/IGVTools/igvtools count -z 5 -w 25 -e 250 "${SAMFILE%%.*}".bam  "${SAMFILE%%.*}".bam.tdf taeGut1

# at first shot, samtools generated error message
# Line 1090, sequence length 1062 vs 1061 from CIGAR
# Parse error at line 1090: CIGAR and sequence length are inconsistent

# look at the offending line
sed -n 1090p 30-gmap-to-tg/lu_master500.sam

####
# trying to fix incorrect CIGAR stings
####

# pysam crashes on encountering incorrect cigar string?
# it's necessary to parse the CIGAR to get expected sequence length
#
# excerpt from samtools source code
#
cat <<CODE
    if (op == 'M') op = BAM_CMATCH;
    else if (op == 'I') op = BAM_CINS;
    else if (op == 'D') op = BAM_CDEL;
    else if (op == 'N') op = BAM_CREF_SKIP;
    else if (op == 'S') op = BAM_CSOFT_CLIP;
    else if (op == 'H') op = BAM_CHARD_CLIP;
    else if (op == 'P') op = BAM_CPAD;
    else if (op == '=') op = BAM_CEQUAL;
    else if (op == 'X') op = BAM_CDIFF; 
    
 int32_t bam_cigar2qlen(const bam1_core_t *c, const uint32_t *cigar)
{
	uint32_t k;
	int32_t l = 0;
	for (k = 0; k < c->n_cigar; ++k) {
		int op = cigar[k] & BAM_CIGAR_MASK;
		if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP || op == BAM_CEQUAL || op == BAM_CDIFF)
			l += cigar[k] >> BAM_CIGAR_SHIFT;
	}
	return l;
} 
CODE

# len = sum(M, I, S, =, X)
# output line indices of lines with CIGAR string not matching the sequence length
gawk -F "\t" '(NF < 10) {print "0 == 0";} (NF > 10) {gsub(/[0-9]+[^MISX=0-9]/, "", $6); gsub(/[MISX=]/, "+", $6); sub(/.$/, " == ", $6); if(length($6) > 4) print $6 length($10); else print "0 == 0"}' 30-gmap-to-tg/lu_master500.sam|bc|grep -n 0 
# to get the raw line numbers of the offending lines
|sed s/:0// > 30-gmap-to-tg/lu_master500.sam.badlines

# first approach - get rid of the bad lines
# (construct sed -e Xd -e Yd ... command line from the number per line input file and run it)
sed 's/^.*$/-e &d/' $SAMFILE.badlines | xargs sed $SAMFILE > ${SAMFILE%%.*}_f.sam
SAMFILE=${SAMFILE%%.*}_f.sam

# the complement to above -- extracting the offending lines
sed 's/^.*$/-e &p/' $SAMFILE.badlines | xargs sed -n $SAMFILE > ${SAMFILE%%.*}_badlines.sam
#------------------------------------------------

# run sim4db with parameters -minXX to report multiple matches to find paralogous genes
# this one is pretty slow -- probably the first process worth parallelization
# or rather prescreening with a fast mapper and using a 'script'
sim4db -genomic /data/genomes/taeGut1/taeGut1.fa -cdna 0a-jp-newbler-contigs/lu_master500.fasta -output 31-sim4db-tg/lu_master500.gff -gff3 -interspecies -mincoverage 70 -minidentity 90 -minlength 60 -alignments -threads 7

#------------------------------------------------
####
# exonerate
####
cat <<QUOTE
I frequently use exonerate to align ESTs to a genome using exonerate server. For example

 fasta2esd genome.fasta genome.esd
 esd2esi genome.esd genome.esi --memorylimit 2048

 exonerate-server genome.esi --port 12886 &

 exonerate --model est2genome ESTs.fasta localhost:12886 --percent 70 --score 100 --showvulgar yes --softmaskquery no --softmasktarget yes --minintron 20 --maxintron 6000 --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --showalignment no --showtargetgff yes --geneseed 250 > $file.exonerate

 Make sure to set percent and maxintron to something sensible for your situation. You can later parse the exonerate output to get GFF. I guess it would also be possible to write a parser to get SAM formatted output but this I don't have.

 The main limitation I find with exonerate is that it's not multi threaded so you have to write your own wrapper to split a large input set of EST and distribute them across threads / nodes.
QUOTE
