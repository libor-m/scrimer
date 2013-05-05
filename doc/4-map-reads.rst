Map reads to reference assembly and detect variants
===================================================

Map contigs to scaffold
-----------------------

Map all the reads using smalt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Set up variables ::

    # data from previous steps
    SCAFFOLD=33-scaffold/lx4.fasta
    INFILES=12-cutadapt/*.fastq

    OUT=40-map-smalt
    SMALT_IDX=${SCAFFOLD%/*}/smalt/k13s4
    CPUS=8

Create index for the scaffold and map the reads::

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

Merge mapping output to a single file 
------------------------------------- 

Create a fasta index for the scaffold::

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

    Copy-paste does not work when copy-pasting to ``bash``, 
    because shell interprets the tab character as autocomplete. Use ``cat`` 
    or escape characters instead.

.. note::

    The order of the rows matters for the vcf output,
    the sample columns order is probably the order of first apperance in the @RG.


Following code generates most of the ``readgroups.txt`` file, you 
have to reorder lines and fill the places marked with '??'::

    DIR=12-cutadapt
    find $DIR -name '*.fastq'|xargs -n1 basename|sed s/.fastq//|gawk '{OFS="\t";print "@RG", "ID:" $0, "SM:??", "LB:" gensub(/\..*$/,"",$0), "PL:LS454", "DS:??";}' > $OUT/readgroups.txt

Prepate the sam files
^^^^^^^^^^^^^^^^^^^^^
Extract the sequence headers, from first ``.sam`` file (the rest of files should be identical)::

    export SAMFILE=$( echo $OUT/*.sam|gawk '{print $1;}' )
    samtools view -S -t $SCAFFOLD.fai -H $SAMFILE > $OUT/sequences.txt
    cat $OUT/sequences.txt $OUT/readgroups.txt > $OUT/sam-header.txt

``samtools merge`` requires sorted alignments, sort them in parallel. This creates ``.bam`` files 
in the output directory::

    parallel "samtools view -but $SCAFFOLD.fai {} | samtools sort - {.}" ::: $OUT/*.sam

Merge it
^^^^^^^^
Merge all the alignments, not removing the duplcates (http://seqanswers.com/forums/showthread.php?t=6543, 
http://seqanswers.com/forums/showthread.php?t=5424). 
# using /[GH]*.bam to avoid generated files (like alldup.bam) in the expansion

    samtools merge -ru -h $OUT/sam-header.txt - $OUT/[GH]*.bam | samtools sort - $OUT/alldup
    samtools index $OUT/alldup.bam


Visualization
-------------
::

    # unmapped read counts
    parallel 'echo $( cut -f2 {}|grep -c "^4$" ) {}' ::: $OUT/*.sam

    # mapping statistics
    samtools idxstats $OUT/alldup.bam|gawk '{map += $3; unmap += $4;} END {print  unmap/map;}'

    # coverage sums for IGV
    igvtools count -z 5 -w 25 -e 250 $OUT/alldup.bam  $OUT/alldup.bam.tdf ${CONTIGS%.*}.genome

