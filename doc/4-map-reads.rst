Map reads to the scaffold
=========================

Map all the reads using smalt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Set up variables:

.. code-block:: bash

    # data from previous steps
    SCAFFOLD=33-scaffold/sc-demo.fasta
    INFILES=10-cutadapt/*.fastq
    OUT=40-map-smalt

    SMALT_IDX=${SCAFFOLD%/*}/smalt/${SCAFFOLD##*/}-k13s4

Create index for the scaffold and map the reads.
Mapping 3 GB of reads (fastq format) takes ~5 hours in 8 threads on Intel Xeon E5620, 0.5 GB memory
per each mapping. 
This step would benefit from parallelization even on multiple machines (not implemented here).

.. code-block:: bash

    # create smalt index
    mkdir -p ${SMALT_IDX%/*}
    smalt index -s 4 $SMALT_IDX $SCAFFOLD

    # map each file, smalt is multithreaded so feed the files sequentially
    mkdir -p $OUT
    parallel -j 1 "smalt map -n $CPUS -p -f sam -o $OUT/{/.}.sam $SMALT_IDX {}" ::: $INFILES

    # Illumina reads can be maped e.g. with bwa
    bwa index $SCAFFOLD
    parallel -j 1 "bwa mem -t $CPUS $SCAFFOLD {} > $OUT/{/.}.sam" ::: $INFILES

Merge mapping output to single file 
-----------------------------------

Create a fasta index for the scaffold:

.. code-block:: bash

    samtools faidx $SCAFFOLD

Create readgroups.txt
^^^^^^^^^^^^^^^^^^^^^

According to your sample wet lab details, create a ``readgroups.txt`` file.
Because ``samtools merge -r`` attaches read group to each alignment (line) in the input 
according to the original filename, the format is ($BASENAME is the fastq file name
without suffix, $SAMPLE is your biological sample, ${BASENAME%%.*} is the dna library name,
all ``tab`` separated)::

    @RG	ID:$BASENAME	SM:$SAMPLE	LB:${BASENAME%%.*}	PL:LS454 DS:$SPECIES

The library name (LB) is important because of ``rmdup``,
description (DS) is here used to identify the species.

.. note::

    The order of the rows matters for the vcf output,
    the sample columns order is probably the order of first appearance in the @RG.

The following code generates most of the ``readgroups.txt`` file, you 
have to reorder lines and fill the places marked with '??':

.. code-block:: bash

    OUT=40-map-smalt
    DIR=10-cutadapt

    find $DIR -name '*.fastq' | xargs -n1 basename | sed s/.fastq// | awk '{OFS="\t";lb=$0;sub(/\..*$/,"",lb);print "@RG", "ID:" $0, "SM:??", "LB:" lb, "PL:LS454", "DS:??";}' > $OUT/readgroups.txt

    # edit the file (ctrl-o enter ctrl-x saves and exits the editor)
    nano $OUT/readgroups.txt

    # sort the readgroups according to species
    <$OUT/readgroups.txt sort -k6,6 > $OUT/readgroups-s.txt

Prepare the sam files
^^^^^^^^^^^^^^^^^^^^^
Extract the sequence headers from the first ``.sam`` file (other files should have identical headers):

.. code-block:: bash

    SAMFILE=$( echo $OUT/*.sam | awk '{print $1;}' )
    samtools view -S -t $SCAFFOLD.fai -H $SAMFILE > $OUT/sequences.txt
    cat $OUT/sequences.txt $OUT/readgroups-s.txt > $OUT/sam-header.txt

``samtools merge`` requires sorted alignments, sort them in parallel. This creates ``.bam`` files 
in the output directory:

.. code-block:: bash

    parallel -j $CPUS "samtools view -but $SCAFFOLD.fai {} | samtools sort - {.}" ::: $OUT/*.sam

Merge
^^^^^
Merge all the alignments. Do not remove duplicates because the duplicate
detection algorithm is based on the read properties of genomic DNA ([#]_, [#]_). 

.. code-block:: bash

    samtools merge -ru -h $OUT/sam-header.txt - $OUT/*.bam | samtools sort - $OUT/alldup
    samtools index $OUT/alldup.bam


Check the results
-----------------

Unmapped read counts.

.. code-block:: bash

    parallel -j $CPUS 'echo $( cut -f2 {}|grep -c "^4$" ) {}' ::: $OUT/*.sam

Mapping statistics

.. code-block:: bash

    samtools idxstats $OUT/alldup.bam | awk '{map += $3; unmap += $4;} END {print  unmap/map;}'

Coverage sums for IGV

.. code-block:: bash

    igvtools count -z 5 -w 25 -e 250 $OUT/alldup.bam  $OUT/alldup.bam.tdf ${CONTIGS%.*}.genome

.. [#] http://seqanswers.com/forums/showthread.php?t=6543 
.. [#] http://seqanswers.com/forums/showthread.php?t=5424
