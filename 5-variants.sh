#------------------------------------------------
# Operation
#------------------------------------------------

# tools used: samtools, tabix, sort-alt, parallel

#  1. call variants
#  2. filter variants

TOOLS=/opt/samtools-0.1.18:/opt/tabix:sort-alt:parallel
export PATH=$TOOLS:$PATH

# 1. call variants
#------------------------------------------------

CONTIGS=33-scaffold/lx4.fasta
ALIGNS=40-map-smalt/alldup.bam
OUT=50-variants
VARIANTS=$OUT/lx4-variants
mkdir -p $OUT
vcfutils.pl splitchr $CONTIGS.fai | parallel "samtools mpileup -DSu -L 10000 -f $CONTIGS -r {} $ALIGNS | bcftools view -bvcg - > $OUT/part-{}.bcf"

# merge the intermediate results
# this generates the correct ordering of the parts, according to the .fai
# max command line size probably should not be a limit since linux 2.6.23 (http://www.in-ulm.de/~mascheck/various/argmax/)
vcfutils.pl splitchr $CONTIGS.fai | xargs -i echo $OUT/part-{}.bcf | xargs -x bcftools cat > $VARIANTS-raw.bcf
bcftools index $VARIANTS-raw.bcf

#TODO: check varFilter options a bit more
bcftools view $VARIANTS-raw.bcf | vcfutils.pl varFilter | bgzip > $VARIANTS-filtered.vcf.gz
tabix -p vcf $VARIANTS-filtered.vcf.gz

# remove the intermediate results, if the merge was ok
vcfutils.pl splitchr $CONTIGS.fai | xargs -i echo $OUT/part-{}.bcf | xargs rm

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


# 2. filter variants
#------------------------------------------------

# filter on average read depth and site quality
# progress meter
VCFINPUT=51-variants-parallel/lx3-variants-filtered.vcf.gz
VCFOUTPUT=51-variants-parallel/lx3-variants-filt2.vcf.gz
pv -p $VCFINPUT | bgzip -d | vcf_filter.py --local-script pyvcf_filters.py --no-filtered - avg-dps sq| bgzip > $VCFOUTPUT
tabix -p vcf $VCFOUTPUT

# filter depth per sample, keep the filtered variants in output
# with progress meter (only adding flags to rows, suppose similar filesize)
VCFINPUT=51-variants-parallel/lx3-variants-filt2.vcf.gz
VCFOUTPUT=51-variants-parallel/lx3-variants-dps.vcf.gz
pv -p $VCFINPUT | bgzip -d | vcf_filter.py --local-script pyvcf_filters.py - dps --depth-per-sample 3 | bgzip > $VCFOUTPUT
tabix -p vcf $VCFOUTPUT

VCFOUTPUT=51-variants-parallel/lx3-variants-selected.vcf.gz
pv -p $VCFINPUT | bgzip -d | vcf_filter.py --local-script pyvcf_filters.py - dps --depth-per-sample 3 snp-only contrast-samples --sample-names ll-lu02 ll-lu05 ll-lu07 | bgzip > $VCFOUTPUT
tabix -p vcf $VCFOUTPUT

# count remaining variants
zcat -d $VCFOUTPUT | grep -c PASS

#------------------------------------------------
# Visualization
#------------------------------------------------

# extract qualities, then check distribution in R (-> common power law distribution, peak at 999)
zcat $VCFINPUT|grep -v '^#'|cut -f6 | R > $VCFINPUT.qual

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
