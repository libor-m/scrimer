#------------------------------------------------
# Operation
#------------------------------------------------

# tools used: samtools, bedtools, tabix, primer3

TOOLS=/opt/samtools-0.1.18:/opt/tabix
export PATH=$TOOLS:$PATH

# data from previous steps
VARIANTS=50-variants/lx4-variants-selected.vcf.gz
CONTIGS=33-scaffold/lx4.fasta
ANNOTS=33-scaffold/lx4.sorted.gff3.gz 

# other inputs/outputs
GFFILE=lx4-primers.gff3
OUT=60-gff-primers
GFF=$OUT/$GFFILE
mkdir -p $OUT

# for all selected variants design pcr and genotyping primers
# takes about a minute for 1000 selected variants, 5 MB gzipped vcf, 26 MB uncompressed genome, 5 MB gzipped gff
./design_primers.py $CONTIGS $ANNOTS $VARIANTS > $GFF

# it's quite vital to sort and index the file to use it in IGV
sortBed -i $GFF | bgzip > ${GFF%.*}.sorted.gff3.gz
tabix -f -p gff ${GFF%.*}.sorted.gff3.gz

# first, check primers with sole blat
# the next step is to set up a gfServer and use gfPcr

#------------------------------------------------
# Visualization
#------------------------------------------------
grep primer-gt $GFF | grep -v -c 'PROBLEMS='
