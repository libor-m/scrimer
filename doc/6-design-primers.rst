.. _primers:

Design primers
==============
Design primers with Primer3
---------------------------
Set inputs and outputs for this step:

.. code-block:: bash

    VARIANTS=50-variants/lx4-variants-selected.vcf.gz
    SCAFFOLD=33-scaffold/lx4.fasta
    ANNOTS=33-scaffold/lx4.sorted.gff3.gz 

    GFFILE=lx4-primers.gff3
    OUT=60-gff-primers
    
    GFF=$OUT/$GFFILE
    PRIMERS=${GFF%.*}.sorted.gff3.gz
    mkdir -p $OUT

    # for all selected variants design pcr and genotyping primers
    # takes about a minute for 1000 selected variants, 5 MB gzipped vcf, 26 MB uncompressed genome, 5 MB gzipped gff
    design_primers.py $SCAFFOLD $ANNOTS $VARIANTS > $GFF

Sort and index the annotation before using it in IGV:

.. code-block:: bash

    sortBed -i $GFF | bgzip > $PRIMERS
    tabix -f -p gff $PRIMERS

Create a region list for IGV to quickly inspect all the primers.

.. code-block:: bash

    awk 'BEGIN{OFS="\t";} /pcr-product/ {print $1, $4, $5;}' $GFF > ${GFF%.*}.bed
    
Convert scaffold to blat format
-------------------------------

.. code-block:: bash

    SCAFFOLD2BIT=$OUT/${SCAFFOLD##*/}.2bit
    faToTwoBit $SCAFFOLD $SCAFFOLD2BIT
    
Validate primers with blat/isPcr
--------------------------------

.. note::
    
    Ideal case would be to add the annotations found by ``isPcr`` to the primer gff3 tags.

Recomended parameters for PCR primers in blat [#]_: ``-tileSize=11``, ``-stepSize=5``

Get the primes sequences, in both formats:
    
.. code-block:: bash

    PRIMERS_BASE=${PRIMERS%%.*}
    extract_primers.py $PRIMERS isPcr > $PRIMERS_BASE.isPcr
    extract_primers.py $PRIMERS > $PRIMERS_BASE.fa

Check against transcriptome data and reference genome:

.. code-block:: bash
    
    # select one of those a time:
    TARGET=$SCAFFOLD2BIT
    TARGET=$GENOME2BIT

    TARGET_TAG=$( basename ${TARGET%%.*} )
    isPcr -out=psl $TARGET $PRIMERS_BASE.isPcr $PRIMERS_BASE.isPcr.$TARGET_TAG.psl
    blat -minScore=15 -tileSize=11 -stepSize=5 -maxIntron=0 $TARGET $PRIMERS_BASE.fa $PRIMERS_BASE.$TARGET_TAG.psl

Check the results
-----------------

See all places where ``primer3`` reported problems:

.. code-block:: bash

    grep primer-gt $GFF | grep -v -c 'PROBLEMS='

Use agrep to find similar sequences in transcript scaffold, to check if the 
settings of blat are ok. Line wrapping in ``fasta`` can lead to false negatives,
but at least some sequences should be found:

.. code-block:: bash

    # agrep is quite enough for simple checks on assemblies of this size (30 MB)
    SEQ=GCACATTTCATGGTCTCCAA
    agrep $SEQ $SCAFFOLD|grep $SEQ

.. [#] http://genomewiki.ucsc.edu/index.php/Blat-FAQ