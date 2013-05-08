Detect and choose variants
==========================

Call variants
-------------
Set up input and output for current step:

.. code-block:: bash

    SCAFFOLD=33-scaffold/lx4.fasta
    ALIGNS=40-map-smalt/alldup.bam

    # outputs
    OUT=50-variants
    VARIANTS=$OUT/lx4-variants
    mkdir -p $OUT

Run the variant calling in parallel. Takes 3 hours for 15 samples on single Intel Xeon E5620:

.. code-block:: bash

    vcfutils.pl splitchr $SCAFFOLD.fai | parallel "samtools mpileup -DSu -L 10000 -f $SCAFFOLD -r {} $ALIGNS | bcftools view -bvcg - > $OUT/part-{}.bcf"

Merge the intermediate results. ``vcfutils.pl`` is used to generate the correct ordering of parts, as in the ``.fai``:

.. code-block:: bash

    vcfutils.pl splitchr $SCAFFOLD.fai | xargs -i echo $OUT/part-{}.bcf | xargs -x bcftools cat > $VARIANTS-raw.bcf
    bcftools index $VARIANTS-raw.bcf

Remove the intermediate results, if the merge was ok:

.. code-block:: bash

    vcfutils.pl splitchr $SCAFFOLD.fai | xargs -i echo $OUT/part-{}.bcf | xargs rm


Filter variants
---------------

We're interested in two kinds of variant qualities 

- all possible variants so they can be avoided in primer design
- high confidence variants that can be used to answer our questions

Filtering strategy:
 
- use common ``samtools`` filtering
  
  - indels caused by 454 homopolymer problems have generally low quality score,
    so they should be filtered at this stage

- remove uninteresting information (for convinient viewing in IGV)
  
  - overall low coverage sites (less than 3 reads per sample - averaged, to avoid discarding
    some otherwise interesting information because of one bad sample)
    
- select the interesting variants, leave the rest in the file flagged as 'uninteresting'
  
  - only SNPs
  - at least 3 reads per sample
  - no shared variants between the two species

Samtools filtering
^^^^^^^^^^^^^^^^^^

We can expect quite high *strand bias* in RNASeq data, so don't filter on strand bias
(``-1 0``), use the defaults for other settings of ``vcfutils varFilter`` command:

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

Use ``pv`` as progress meter. ``pv`` can be substituted by ``cat``:

.. code-block:: bash

    # filter on average read depth and site quality
    VCFINPUT=$VARIANTS-filtered.vcf.gz
    VCFOUTPUT=$VARIANTS-filt2.vcf.gz
    pv -p $VCFINPUT | bgzip -d | vcf_filter.py --no-filtered - avg-dps sq| bgzip > $VCFOUTPUT
    tabix -p vcf $VCFOUTPUT

Interesting variants
^^^^^^^^^^^^^^^^^^^^

Keep the rest in the file, with mark in ``FILTER`` filed. This way both 
the selected and non-selected variants can be checked in IGV:

.. code-block:: bash

    VCFINPUT=$VARIANTS-filt2.vcf.gz
    VCFOUTPUT=$VARIANTS-selected.vcf.gz
    pv -p $VCFINPUT | bgzip -d | vcf_filter.py - dps --depth-per-sample 3 snp-only contrast-samples --sample-names lu02 lu05 lu07 lu10 lu12 lu14 lu15| bgzip > $VCFOUTPUT
    tabix -p vcf $VCFOUTPUT

Check the results
-----------------

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

