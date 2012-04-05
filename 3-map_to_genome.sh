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

# build gmap index for zebra finch
gmap_build -d gmap_taeGut1 -D /data/genomes/taeGut1 /data/genomes/taeGut1/taeGut1.fa

# map the jpaces contigs to the zebra finch genome
gmap -D /data/genomes/taeGut1 -d gmap_taeGut1 -f samse -B 3 -x 30 -t 6 --cross-species --read-group-id=lu_master500 0a-jp-newbler-contigs/lu_master500.fasta  > 30-gmap-to-tg/lu_master500.sam

SAMFILE=30-gmap-to-tg/lu_master500.sam

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

# map again, using gff3_gene as output for interoperability with sim4db
# takes 25 mins, better than trying to convert sam to gff
gmap -D /data/genomes/taeGut1 -d gmap_taeGut1 -f gff3_gene -B 3 -x 30 -t 8 --cross-species 0a-jp-newbler-contigs/lu_master500.fasta  > 30-gmap-to-tg/lu_master500.gff

####
# sim4db
####
# run sim4db with parameters -minXX to report multiple matches to find paralogous genes
# this one is pretty slow -- probably the first process worth parallelization
# or rather prescreening with a fast mapper and using a 'script'
sim4db -genomic /data/genomes/taeGut1/taeGut1.fa -cdna 0a-jp-newbler-contigs/lu_master500.fasta -output 31-sim4db-tg/lu_master500.gff -gff3 -interspecies -mincoverage 70 -minidentity 90 -minlength 60 -alignments -threads 7

# use a fast mapper to find all +-50 KB windows for predicting exon/gene models with sim4db
# smalt with -d on shattered contigs could work reasonably well and fast

# create fragments, using slightly modified fasta_fragments.py from lastz distribution
cat 0a-jp-newbler-contigs/lu_master500.fasta | ./fasta_fragments.py --step=80 > 31-sim4db-tg/lu_master500_frags.fasta

# map the fragments with smalt (takes only 6 minutes!), reporting all hits (-d -1) scoring over 60
smalt map -n 8 -f cigar -o 31-sim4db-tg/lu_master500_frags.cigar -d -1 -m 60 /data/genomes/taeGut1/smalt/taeGut1k13s4 31-sim4db-tg/lu_master500_frags.fasta

# construct the script for sim4db
# this could be probably achieved also by pipeline of bedtools (slopBed|mergeBed|sed/awk) -- should be done for publication?
# (maybe the python script is more comprehensible)
cat 31-sim4db-tg/lu_master500_frags.cigar|./cigar_to_sim4db_scr.py /data/genomes/taeGut1/taeGut1.fa.fai | sort --key=5 -n > 31-sim4db-tg/lu_master500.sim4scr

# run sim4db using the script (takes several seconds for the whole genome !!)
sim4db -genomic /data/genomes/taeGut1/taeGut1.fa -cdna 0a-jp-newbler-contigs/lu_master500.fasta -script 31-sim4db-tg/lu_master500.sim4scr -output 31-sim4db-tg/lu_master500_scr.gff -gff3 -interspecies -mincoverage 70 -minidentity 90 -minlength 60 -alignments -threads 7

# extract chromosome 10 for comparison with the non-seeded method
grep chr10 31-sim4db-tg/lu_master500_scr.gff|sed s/0:chr/chr/ > 31-sim4db-tg/test_chr10_scr.gff

# compare the annotations 
# to be done programaticaly

# fix chromosome names for the whole gff
cat 31-sim4db-tg/lu_master500_scr.gff|sed s/^[0-9][0-9]*:chr/chr/ > 31-sim4db-tg/lu_master500_scr_fix.gff

# run the fixed sim4db again
sim4db -genomic /data/genomes/taeGut1/taeGut1.fa -cdna 0a-jp-newbler-contigs/lu_master500.fasta -script 31-sim4db-tg/lu_master500.sim4scr -output 31-sim4db-tg/lu_master500_sim4db2.gff -gff3 -interspecies -mincoverage 70 -minidentity 90 -minlength 60 -alignments -threads 7

# fix chromosome names 
sed s/^[0-9][0-9]*:chr/chr/ 31-sim4db-tg/lu_master500_sim4db2.gff > 31-sim4db-tg/lu_master500_sim4db2_fix.gff

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

####
# transfer the annotation
####

# select features in EnsGene file 
# probably a dead end, the information will be transferred during liftover
intersectBed -wa -a /data/genomes/taeGut1/annot/ensGene_s.bed.gz -b 31-sim4db-tg/lu_master500_scr_fix.gff | uniq > 32-extracted-annotation/ensGene_sim4.bed
intersectBed -wa -a /data/genomes/taeGut1/annot/ensGene_s.bed.gz -b 30-gmap-to-tg/lu_master500.gff | uniq > 32-extracted-annotation/ensGene_gmap.bed
cat 32-extracted-annotation/ensGene_*.bed|sortBed -i stdin|uniq > 32-extracted-annotation/ensGene.bed

intersectBed -wa -a /data/genomes/taeGut1/annot/xenoMrna_s.bed.gz -b 30-gmap-to-tg/lu_master500.gff | uniq > 32-extracted-annotation/xenoMrna_gmap.bed
intersectBed -wa -a /data/genomes/taeGut1/annot/xenoMrna_s.bed.gz -b 31-sim4db-tg/lu_master500_scr_fix.gff | uniq > 32-extracted-annotation/xenoMrna_sim4.bed

# liftover using the exon mapper coordinates
#TODO: one of the mappers inverts the 'Target' coordinates when target and strand are both -
# have a look at gff spec and find the solution
#
# sim4db manual (http://sourceforge.net/apps/mediawiki/kmer/index.php?title=Getting_Started_with_Sim4db)
# Exon coordinates are nucleotide based, starting from 1. Genomic coordinates are always in the original sequence, 
# while the cDNA coordinates will refer to positions in the reverse complement of the sequence if the match orientation is indicated as 'complement'.
#
# --> this is unnecessary, because the orientation of the transcript can be deduced from the target chromosome strand ..?
# --> patched in sim4db, patch on sourceforge

# use liftover to transfer the ensGenes exon annotations
./liftover.py 31-sim4db-tg/lu_master500_sim4db2_fix.gff /data/genomes/taeGut1/annot/ensGene_s.bed.gz > 32-liftover/sim4db-ensGenes.gff
./liftover.py 30-gmap-to-tg/lu_master500.gff /data/genomes/taeGut1/annot/ensGene_s.bed.gz > 32-liftover/gmap-ensGenes.gff

# construct a virtual genome (contigs joined in order of appearance on reference genome chromosomes)
# 'N' gaps should be larger than max read size to avoid the mapping of the reads across gaps
./virtual_genome.py 0a-jp-newbler-contigs/lu_master500.fasta 32-liftover/gmap-ensGenes.gff 32-liftover/sim4db-ensGenes.gff 33-virtual-genome/lx2.fasta 33-virtual-genome/lx2.gff
./virtual_genome.py 0a-jp-newbler-contigs/lu_master500.fasta 32-liftover/gmap-ensGenes.gff 32-liftover/sim4db-ensGenes.gff 33-virtual-genome/lx3.fasta 33-virtual-genome/lx3.gff3

# it's quite vital to sort and index the file to use it in IGV
sortBed -i 33-virtual-genome/lx3.gff3 > 33-virtual-genome/lx3.sorted.gff3
bgzip 33-virtual-genome/lx3.sorted.gff3
tabix -p gff 33-virtual-genome/lx3.sorted.gff3.gz
