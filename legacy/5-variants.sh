#------------------------------------------------
# Operation
#------------------------------------------------

# tools used: samtools, tabix, parallel, pv

#  1. call variants
#  2. filter variants

TOOLS=/opt/samtools-0.1.18:/opt/tabix:/opt/bedtools/bin
export PATH=$TOOLS:$PATH

# 1. call variants
#------------------------------------------------

# data from previous steps
CONTIGS=33-scaffold/lx4.fasta
ALIGNS=40-map-smalt/alldup.bam

# outputs
OUT=50-variants
VARIANTS=$OUT/lx4-variants
mkdir -p $OUT

# takes 3 hours for 15 samples on single Intel Xeon E5620
vcfutils.pl splitchr $CONTIGS.fai | parallel "samtools mpileup -DSu -L 10000 -f $CONTIGS -r {} $ALIGNS | bcftools view -bvcg - > $OUT/part-{}.bcf"

# merge the intermediate results
# this generates the correct ordering of the parts, according to the .fai
# max command line size probably should not be a limit since linux 2.6.23 (http://www.in-ulm.de/~mascheck/various/argmax/)
vcfutils.pl splitchr $CONTIGS.fai | xargs -i echo $OUT/part-{}.bcf | xargs -x bcftools cat > $VARIANTS-raw.bcf
bcftools index $VARIANTS-raw.bcf

# varFilter -1 is strand bias, we can expect quite high strand bias in 
# rnaseq data..
bcftools view $VARIANTS-raw.bcf | vcfutils.pl varFilter -1 0 | bgzip > $VARIANTS-filtered.vcf.gz
tabix -p vcf $VARIANTS-filtered.vcf.gz

# remove the intermediate results, if the merge was ok
vcfutils.pl splitchr $CONTIGS.fai | xargs -i echo $OUT/part-{}.bcf | xargs rm


# 2. filter variants
#------------------------------------------------

# filter the variants
# first completely remove the really uninteresting information
# (for convinient viewing in IGV)
# - overall low coverage sites (less than 3 reads per sample - averaged, to avoid discarding
#   an otherwise interesting information because of one bad sample)
# - indels caused by 454 homopolymer problems 
#   can be recognized by strand bias in DP4 INFO tag
#   INFO=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
#   single base indels with strand bias
#   ! but generally those should have low quality score from samtools, so just filter the low 
#   ! quality variants


# filter on average read depth and site quality
# pv as progress meter
VCFINPUT=$VARIANTS-filtered.vcf.gz
VCFOUTPUT=$VARIANTS-filt2.vcf.gz
pv -p $VCFINPUT | bgzip -d | vcf_filter.py --local-script pyvcf_filters.py --no-filtered - avg-dps sq| bgzip > $VCFOUTPUT
tabix -p vcf $VCFOUTPUT

# filter depth per sample, keep the filtered variants in output
# --sample-names is the SM: tag in readgroups.txt
VCFINPUT=$VARIANTS-filt2.vcf.gz
VCFOUTPUT=$VARIANTS-selected.vcf.gz
pv -p $VCFINPUT | bgzip -d | vcf_filter.py --local-script pyvcf_filters.py - dps --depth-per-sample 3 snp-only contrast-samples --sample-names lu02 lu05 lu07 lu10 lu12 lu14 lu15| bgzip > $VCFOUTPUT
tabix -p vcf $VCFOUTPUT

#------------------------------------------------
# Visualization
#------------------------------------------------

# extract qualities, then check distribution in R (-> common power law distribution, peak at 999)
zcat $VCFINPUT|grep -v '^#'|cut -f6 | R > $VCFINPUT.qual

# count remaining variants
zcat -d $VCFOUTPUT | grep -c PASS

# count variants on chromosome Z
zcat -d $VCFOUTPUT | grep PASS | grep -c ^chrZ


#------------------------------------------------
# spare parts
#------------------------------------------------

# single file workflow
#------------------------------------------------
# do the genotype likelihood computations
# takes about three hours on Intel Xeon E5620, produces about 4 GB of data (I shouldn't have used the -u)
samtools mpileup -D -S -u -L 10000 -f 33-virtual-genome/lx3.fasta 41-vg-map-smalt/alldup.bam > 50-variants/lx3-mpileup.bcf

# call the genotypes
bcftools view -bvcg 50-variants/lx3-mpileup.bcf > 50-variants/lx3-variants-raw.bcf  

# filter and output vcf
bcftools view 50-variants/lx3-variants-raw.bcf | vcfutils.pl varFilter > 50-variants/lx3-variants-filtered.vcf

# index with tabix to use in igv
bgzip 50-variants/lx3-variants-filtered.vcf
tabix -p vcf 50-variants/lx3-variants-filtered.vcf.gz

# filter 
VCFINPUT=$VARIANTS-filt2.vcf.gz
VCFOUTPUT=$VARIANTS-dps.vcf.gz
pv -p $VCFINPUT | bgzip -d | vcf_filter.py --local-script pyvcf_filters.py - dps --depth-per-sample 3 | bgzip > $VCFOUTPUT
tabix -p vcf $VCFOUTPUT
