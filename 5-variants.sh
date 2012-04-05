TOOLS=samtools:tabix:sort-alt:parallel

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

#
# parallel workflow (less than one hour)
#
CONTIGS=33-virtual-genome/lx3.fasta
ALIGNS=41-vg-map-smalt/alldup2.bam
OUTDIR=51-variants-parallel
vcfutils.pl splitchr $CONTIGS.fai | parallel "samtools mpileup -DSu -L 10000 -f $CONTIGS -r {} $ALIGNS | bcftools view -bvcg - > $OUTDIR/part-{}.bcf"

# merge the intermediate results
# this generates the correct ordering of the parts, according to the .fai
# max command line size probably should not be a limit since linux 2.6.23 (http://www.in-ulm.de/~mascheck/various/argmax/)
vcfutils.pl splitchr $CONTIGS.fai | xargs -i echo $OUTDIR/part-{}.bcf | xargs -x bcftools cat > $OUTDIR/lx3-variants-raw.bcf
bcftools index $OUTDIR/lx3-variants-raw.bcf

bcftools view $OUTDIR/lx3-variants-raw.bcf | vcfutils.pl varFilter | bgzip > $OUTDIR/lx3-variants-filtered.vcf.gz
tabix -p vcf $OUTDIR/lx3-variants-filtered.vcf.gz

# remove the intermediate results, if the merge was ok
vcfutils.pl splitchr $CONTIGS.fai | xargs -i echo $OUTDIR/part-{}.bcf | xargs rm

# filter the variants
# first completely remove the really uninteresting information
# (for convinient viewing in IGV)
# - overall low coverage sites (less than 3 reads per sample - averaged, to avoid discarding
#   an otherwise interesting information because of one bad sample)
# - indels caused by 454 homopolymer problems - can be recognized by strand bias in DP4 INFO tag
#   INFO=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
#   single base indels with strand bias

# extract qualities, then check distribution in R (-> common power law distribution, peak at 999)
zcat $VCFINPUT|grep -v '^#'|cut -f6 | R > $VCFINPUT.qual

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
vcf_filter.py --local-script pyvcf_filters.py --depth-per-sample 3 $VCFINPUT dps | bgzip | pv -s $( stat -c%s $VCFINPUT ) > $VCFOUTPUT
tabix -p vcf $VCFOUTPUT

