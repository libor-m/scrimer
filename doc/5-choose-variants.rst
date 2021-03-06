Detect and choose variants
==========================

Call variants with samtools
---------------------------

Set up input and output for the current step:

.. code-block:: bash

    SCAFFOLD=33-scaffold/sc-demo.fasta
    ALIGNS=40-map-smalt/alldup.bam

    # outputs
    OUT=50-variants
    VARIANTS=$OUT/demo-variants
    mkdir -p $OUT

Run the variant calling in parallel. Takes 3 hours for 15 samples on a single Intel Xeon E5620:

.. code-block:: bash

    vcfutils.pl splitchr $SCAFFOLD.fai | parallel -j $CPUS "samtools mpileup -DSu -L 10000 -f $SCAFFOLD -r {} $ALIGNS | bcftools view -bvcg - > $OUT/part-{}.bcf"

    # samtools > 0.1.19
    vcfutils.pl splitchr $SCAFFOLD.fai | parallel -j $CPUS "samtools mpileup -u -t DP,SP -L 10000 -f $SCAFFOLD -r {} $ALIGNS | bcftools call -O b -vm - > $OUT/part-{}.bcf"

Merge the intermediate results. ``vcfutils.pl`` is used to generate the correct ordering of parts, as in the ``.fai``:

.. code-block:: bash

    vcfutils.pl splitchr $SCAFFOLD.fai | xargs -i echo $OUT/part-{}.bcf | xargs -x bcftools cat > $VARIANTS-raw.bcf
    bcftools index $VARIANTS-raw.bcf

    # samtools > 0.1.19
    # workaround for https://github.com/samtools/bcftools/issues/134
    vcfutils.pl splitchr $SCAFFOLD.fai | parallel "bcftools view -O z $OUT/part-{}.bcf > $OUT/{.}.vcf.gz && bcftools index $OUT/{.}.vcf.gz"
    bcftools concat -O b $OUT/*.vcf.gz > $VARIANTS-raw.bcf
    bcftools index $VARIANTS-raw.bcf

Remove the intermediate results, if the merge was ok:

.. code-block:: bash

    vcfutils.pl splitchr $SCAFFOLD.fai | xargs -i echo $OUT/part-{}.bcf | xargs rm

Call variants with FreeBayes
----------------------------

This is an alternative to the previous section. FreeBayes uses local realignment around INDELs, so the 
variant calling for 454 data should be better.

.. code-block:: bash

    SCAFFOLD=33-scaffold/sc-demo.fasta
    ALIGNS=40-map-smalt/alldup.bam
    OUT=51-var-freebayes

    GNAME=$( echo ${SCAFFOLD##*/} | cut -d. -f1 )
    mkdir -p $OUT
    vcfutils.pl splitchr $SCAFFOLD.fai | parallel -j $CPUS "freebayes -f $SCAFFOLD -r {} $ALIGNS > $OUT/${GNAME}-{}.vcf"

    # join the results
    OFILE=$OUT/variants-raw.vcf

    # headers
    FILE=$( find $OUT -name ${GNAME}-*.vcf | sort | head -1 )
    <$FILE egrep '^#' > $OFILE
    # the rest in .fai order
    vcfutils.pl splitchr $SCAFFOLD.fai | parallel -j 1 "egrep -v '^#' $OUT/${GNAME}-{}.vcf >> $OFILE"
    
    # filter the variants on quality (ignore the warning messages)
    <$OFILE bcftools view -i 'QUAL > 20' -O z > $OUT/variants-qual.vcf.gz
    tabix -p vcf $OUT/variants-qual.vcf.gz

Filter variants
---------------

We're interested in two kinds of variant qualities 

- all possible variants so they can be avoided in primer design
- high confidence variants that can be used to answer our questions

Filtering strategy:
 
- use predefined ``samtools`` filtering
  
  - indels caused by 454 homopolymer problems generally have low quality scores,
    so they should be filtered at this stage

- remove uninteresting information (for convenient viewing in IGV)
  
  - overall low coverage sites (less than 3 reads per sample - averaged, to avoid discarding
    some otherwise interesting information because of one bad sample)
    
- select the interesting variants, leave the rest in the file flagged as 'uninteresting'
  
  - only SNPs
  - at least 3 reads per sample
  - no shared variants between the two species

Samtools filtering
^^^^^^^^^^^^^^^^^^

Quite high *strand bias* in RNASeq data can be expected, so don't filter on strand bias
(``-1 0``). Use the defaults for other settings of ``vcfutils varFilter`` command:

- minimum RMS mapping quality for SNPs [10]
- minimum read depth [2]
- maximum read depth [10000000]
- minimum number of alternate bases [2]
- SNP within INT bp around a gap to be filtered [3]
- window size for filtering adjacent gaps [10]
- min P-value for baseQ bias [1e-100]
- min P-value for mapQ bias [0]
- min P-value for end distance bias [0.0001]
- FLOAT  min P-value for HWE (plus F<0) [0.0001]

.. code-block:: bash

    bcftools view $VARIANTS-raw.bcf | vcfutils.pl varFilter -1 0 | bgzip > $VARIANTS-filtered.vcf.gz
    tabix -p vcf $VARIANTS-filtered.vcf.gz

Convenience filtering
^^^^^^^^^^^^^^^^^^^^^
The required average depth per sample can be adjusted here.
Using ``pv`` as a progress meter. ``pv`` can be substituted by ``cat``:

.. code-block:: bash

    # filter on average read depth and site quality
    VCFINPUT=$VARIANTS-filtered.vcf.gz
    VCFOUTPUT=$VARIANTS-filt2.vcf.gz
    pv -p $VCFINPUT | bgzip -d | vcf_filter.py --no-filtered - avg-dps --avg-depth-per-sample 5 sq --site-quality 30 | bgzip > $VCFOUTPUT

    # samtools > 0.1.19 produce conflicting info tags, get rid of it if the above filtering fails
    pv -p $VCFINPUT | bgzip -d | sed 's/,Version="3"//' | vcf_filter.py --no-filtered - avg-dps --avg-depth-per-sample 5 sq --site-quality 30 | bgzip > $VCFOUTPUT

    # freebayes output
    zcat $OUT/variants-qual.vcf.gz| vcfutils.pl varFilter -1 0 | vcf_filter.py --no-filtered - avg-dps --avg-depth-per-sample 5 sq --site-quality 30 | bgzip > $OUT/demo-filt1.vcf.gz

Interesting variants
^^^^^^^^^^^^^^^^^^^^
The filtered out variants are kept in the file, only marked as filtered out. This way both 
the selected and non-selected variants can be checked in IGV. Required minumum depth per sample can be adjusted here:

.. code-block:: bash

    # samtools files
    VCFINPUT=$VARIANTS-filt2.vcf.gz
    VCFOUTPUT=$VARIANTS-selected.vcf.gz
    
    # freebayes files
    VCFINPUT=$OUT/demo-filt1.vcf.gz
    VCFOUTPUT=$OUT/demo-selected.vcf.gz
    
    <$VCFINPUT pv -p | zcat | vcf_filter.py - dps --depth-per-sample 3 snp-only contrast-samples --sample-names lu02 lu14 lu15 | bgzip > $VCFOUTPUT
    tabix -p vcf $VCFOUTPUT
    
Check the results
-----------------

Extract calculated variant  qualities, so the distribution
can be checked (-> common power law distribution, additional peak at 999):

.. code-block:: bash

    zcat $VCFINPUT | grep -v '^#' | cut -f6 > $VCFINPUT.qual
    # and visualize externally

    # or directly in terminal, using bit.ly data_hacks tools
    zcat $VCFINPUT | grep -v '^#' | cut -f6 | histogram.py -b 30

Count selected variants:

.. code-block:: bash

    zcat -d $VCFOUTPUT | grep -c PASS

Count variants on **chromosome Z**:

.. code-block:: bash

    zcat -d $VCFOUTPUT | grep PASS | grep -c ^chrZ

