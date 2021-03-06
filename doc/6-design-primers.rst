.. _primers:

Design primers
==============
Design primers with Primer3
---------------------------
Set inputs and outputs for this step:

.. code-block:: bash

    VARIANTS=50-variants/demo-variants-selected.vcf.gz
    SCAFFOLD=33-scaffold/sc-demo.fasta
    ANNOTS=33-scaffold/sc-demo.sorted.gff3.gz

    GFFILE=demo-primers.gff3
    OUT=60-gff-primers
    
    GFF=$OUT/$GFFILE
    PRIMERS=${GFF%.*}.sorted.gff3.gz
    mkdir -p $OUT

    # for all selected variants design pcr and genotyping primers
    # takes about a minute for 1000 selected variants, 5 MB gzipped vcf, 26 MB uncompressed genome, 5 MB gzipped gff
    # default values are set for SNaPshot
    export PRIMER3_CONFIG=/opt/primer3/primer3_config/
    design_primers.py $SCAFFOLD $ANNOTS $VARIANTS > $GFF
    
    # use --primer-pref to set preferred length of genotyping primer
    # this is useful for other genotyping methods, like MALDI-TOF
    design_primers.py --primer-pref 15 --primer-max 25 $SCAFFOLD $ANNOTS $VARIANTS > $GFF

Sort and index the annotation before using it in IGV. For a small set of primers it is not necessary to 
compress and index the file, IGV can handle raw files as well.

.. code-block:: bash

    sortBed -i $GFF | bgzip > $PRIMERS
    tabix -f -p gff $PRIMERS

Create a region list for IGV to quickly inspect all the primers. 

.. code-block:: bash

    <$GFF awk 'BEGIN{OFS="\t";} /pcr-product/ {match($9, "ID=[^;]+"); print $1, $4, $5, substr($9, RSTART+3, RLENGTH);}' > ${GFF%.*}.bed
    
Convert scaffold to the blat format
-----------------------------------

.. code-block:: bash

    SCAFFOLD2BIT=$OUT/${SCAFFOLD##*/}.2bit
    faToTwoBit $SCAFFOLD $SCAFFOLD2BIT
    
Validate primers with blat/isPcr
--------------------------------

Recommended parameters for PCR primers in blat [#]_: ``-tileSize=11``, ``-stepSize=5``

Get the primer sequences, in formats for isPcr and blat:
    
.. code-block:: bash

    PRIMERS_BASE=${PRIMERS%%.*}
    extract_primers.py $PRIMERS isPcr > $PRIMERS_BASE.isPcr
    extract_primers.py $PRIMERS > $PRIMERS_BASE.fa

Check against transcriptome data and the reference genome:

.. code-block:: bash
    
    # select one of those a time:
    TARGET=$SCAFFOLD2BIT
    TARGET=$GENOME2BIT

    TARGET_TAG=$( basename ${TARGET%%.*} )
    isPcr -out=psl $TARGET $PRIMERS_BASE.isPcr $PRIMERS_BASE.isPcr.$TARGET_TAG.psl
    blat -minScore=15 -tileSize=11 -stepSize=5 -maxIntron=0 $TARGET $PRIMERS_BASE.fa $PRIMERS_BASE.$TARGET_TAG.psl

.. note::
    
    TODO: It would be nice to add the annotations found by ``isPcr`` to the primer gff3 tags (not implemented yet). 

Check the results
-----------------

Count and check all places where ``primer3`` reported problems:

.. code-block:: bash

    <$GFF grep gt-primer | grep -c 'PROBLEMS='
    <$GFF grep gt-primer | grep 'PROBLEMS=' | less -S

    # count unique variants with available primer set
    <$GFF grep gt-primer|grep -v PROBLEM|egrep -o 'ID=[^;]+'|cut -c-13|sort -u|wc -l

Use agrep to find similar sequences in the transcript scaffold, to check if the 
sensitivity settings of blat are OK. Line wrapping in ``fasta`` can lead to false negatives,
but at least some primers should yield hits:

.. code-block:: bash

    # agrep is quite enough for simple checks on assemblies of this size (30 MB)
    SEQ=GCACATTTCATGGTCTCCAA
    agrep $SEQ $SCAFFOLD|grep $SEQ

Import your primers to any spreadsheet program with some selected information on each
primer. Use copy and paste, the file format is ``tab`` separated values. When there is more 
than one genotyping primer for one PCR product, the information on the PCR product is repeated.

.. code-block:: bash

    extract_primers.py $PRIMERS table > $PRIMERS_BASE.tsv

.. [#] http://genomewiki.ucsc.edu/index.php/Blat-FAQ