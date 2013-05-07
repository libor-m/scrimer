Choose interesting variants
===========================

Call variants
-------------
Set up input output for current step:

.. code-block:: bash

    CONTIGS=33-scaffold/lx4.fasta
    ALIGNS=40-map-smalt/alldup.bam

    # outputs
    OUT=50-variants
    VARIANTS=$OUT/lx4-variants
    mkdir -p $OUT

Run the variant calling in parallel. Takes 3 hours for 15 samples on single Intel Xeon E5620:

.. code-block:: bash

    vcfutils.pl splitchr $CONTIGS.fai | parallel "samtools mpileup -DSu -L 10000 -f $CONTIGS -r {} $ALIGNS | bcftools view -bvcg - > $OUT/part-{}.bcf"

Merge the intermediate results. ``vcfutils.pl`` is used to generate the correct ordering of parts, as in the ``.fai``:

.. code-block:: bash

    vcfutils.pl splitchr $CONTIGS.fai | xargs -i echo $OUT/part-{}.bcf | xargs -x bcftools cat > $VARIANTS-raw.bcf
    bcftools index $VARIANTS-raw.bcf

``varFilter``'s  ``-1`` argument is strand bias: we can expect quite a high strand bias in 
rnaseq data, so don't filter on strand bias:

.. code-block:: bash

    bcftools view $VARIANTS-raw.bcf | vcfutils.pl varFilter -1 0 | bgzip > $VARIANTS-filtered.vcf.gz
    tabix -p vcf $VARIANTS-filtered.vcf.gz

Remove the intermediate results, if the merge was ok:

.. code-block:: bash

    vcfutils.pl splitchr $CONTIGS.fai | xargs -i echo $OUT/part-{}.bcf | xargs rm


Filter variants
---------------

- completely remove the really uninteresting information (for convinient viewing in IGV)
  
  - overall low coverage sites (less than 3 reads per sample - averaged, to avoid discarding
    some otherwise interesting information because of one bad sample)
  - indels caused by 454 homopolymer problems have generally low quality score,
    so it should not be necessary to filter them separately

Basic filtering. Usew ``pv`` as progress meter. ``pv`` can be substituted by ``cat``:

.. code-block:: bash

    # filter on average read depth and site quality
    VCFINPUT=$VARIANTS-filtered.vcf.gz
    VCFOUTPUT=$VARIANTS-filt2.vcf.gz
    pv -p $VCFINPUT | bgzip -d | vcf_filter.py --local-script pyvcf_filters.py --no-filtered - avg-dps sq| bgzip > $VCFOUTPUT
    tabix -p vcf $VCFOUTPUT

Select the **relevant variants**, keep the rest in the file. This way both 
the selected and non-selected variants can be seen in IGV:

.. code-block:: bash

    VCFINPUT=$VARIANTS-filt2.vcf.gz
    VCFOUTPUT=$VARIANTS-selected.vcf.gz
    pv -p $VCFINPUT | bgzip -d | vcf_filter.py --local-script pyvcf_filters.py - dps --depth-per-sample 3 snp-only contrast-samples --sample-names lu02 lu05 lu07 lu10 lu12 lu14 lu15| bgzip > $VCFOUTPUT
    tabix -p vcf $VCFOUTPUT

Visualization
-------------

Extract calculated variant  qualities, so the distribution
can be checked (-> common power law distribution, additional peak at 999):

.. code-block:: bash

    zcat $VCFINPUT | grep -v '^#' | cut -f6 > $VCFINPUT.qual

Count selected variants:

.. code-block:: bash

    zcat -d $VCFOUTPUT | grep -c PASS

Count variants on **chromosome Z**:

.. code-block:: bash

    zcat -d $VCFOUTPUT | grep PASS | grep -c ^chrZ

