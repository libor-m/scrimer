Map reads to reference assembly
===============================

Map contigs to scaffold
-----------------------

Map all the reads using smalt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Set up variables:

.. code-block:: bash

    # data from previous steps
    SCAFFOLD=33-scaffold/lx4.fasta
    INFILES=12-cutadapt/*.fastq

    OUT=40-map-smalt
    SMALT_IDX=${SCAFFOLD%/*}/smalt/${SCAFFOLD##*/}-k13s4
    CPUS=8

Create index for the scaffold and map the reads.
Mapping 3 GB of reads (fastq format) takes ~5 hours in 8 threads on Intel Xeon E5620, 0.5 GB memory
per each mapping. 
This step is probably worth some parallelization on multiple machines.

.. note::

    we used smalt-0.7.0.1, because smalt-0.7.4 was crashing with 
    
    $ [0] rmap.c:1448 ERROR: assertion failed

.. code-block:: bash

    # create smalt index
    mkdir -p ${SMALT_IDX%/*}
    smalt index -s 4 $SMALT_IDX $SCAFFOLD

    # map each file, smalt is multithreaded so feed the files sequentially
    mkdir -p $OUT
    for FQFILE in $INFILES
    do
      SAMFILE=$OUT/$( basename ${FQFILE%.*} ).sam
      smalt map -n $CPUS -p -f sam -o $SAMFILE $SMALT_IDX $FQFILE
    done

Merge mapping output to single file 
-----------------------------------

Create a fasta index for the scaffold:

.. code-block:: bash

    samtools faidx $SCAFFOLD

Create readgroups.txt
^^^^^^^^^^^^^^^^^^^^^

According to real sample info, create a ``readgroups.txt`` file.
Because ``samtools merge -r`` attaches read group to each alignment (line) in input 
according to original filename, the format is (\t separated)::

    @RG	ID:$BASENAME	SM:$SAMPLE	LB:${BASENAME%%.*}	PL:LS454

The library name (LB) is important because of ``rmdup``,
description (DS) is used to identify the species

.. note::

    The order of the rows matters for the vcf output,
    the sample columns order is probably the order of first apperance in the @RG.

Following code generates most of the ``readgroups.txt`` file, you 
have to reorder lines and fill the places marked with '??':

.. code-block:: bash

    OUT=40-map-smalt
    DIR=12-cutadapt
    find $DIR -name '*.fastq'|xargs -n1 basename|sed s/.fastq//|gawk '{OFS="\t";print "@RG", "ID:" $0, "SM:??", "LB:" gensub(/\..*$/,"",$0), "PL:LS454", "DS:??";}' > $OUT/readgroups.txt

Prepate the sam files
^^^^^^^^^^^^^^^^^^^^^
Extract the sequence headers, from first ``.sam`` file (the rest of files should be identical):

.. code-block:: bash

    export SAMFILE=$( echo $OUT/*.sam|gawk '{print $1;}' )
    samtools view -S -t $SCAFFOLD.fai -H $SAMFILE > $OUT/sequences.txt
    cat $OUT/sequences.txt $OUT/readgroups.txt > $OUT/sam-header.txt

``samtools merge`` requires sorted alignments, sort them in parallel. This creates ``.bam`` files 
in the output directory:

.. code-block:: bash

    parallel "samtools view -but $SCAFFOLD.fai {} | samtools sort - {.}" ::: $OUT/*.sam

Merge it
^^^^^^^^
Merge all the alignments. Do not remove duplcates because the duplicate
detection algorithm is based on read properties of genomic DNA ([#]_, [#]_). 

``/[GH]*.bam`` avoids generated files like ``alldup.bam`` in glob expansion.

.. code-block:: bash

    samtools merge -ru -h $OUT/sam-header.txt - $OUT/[GH]*.bam | samtools sort - $OUT/alldup
    samtools index $OUT/alldup.bam


Check the results
-----------------

Unmapped read counts.

.. code-block:: bash

    parallel 'echo $( cut -f2 {}|grep -c "^4$" ) {}' ::: $OUT/*.sam

Mapping statistics

.. code-block:: bash

    samtools idxstats $OUT/alldup.bam|gawk '{map += $3; unmap += $4;} END {print  unmap/map;}'

Coverage sums for IGV

.. code-block:: bash

    igvtools count -z 5 -w 25 -e 250 $OUT/alldup.bam  $OUT/alldup.bam.tdf ${CONTIGS%.*}.genome

.. [#] http://seqanswers.com/forums/showthread.php?t=6543 
.. [#] http://seqanswers.com/forums/showthread.php?t=5424
