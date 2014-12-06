Prepare reference genome
========================

Downlaod and prepare the reference genome
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- List of available genomes is at http://hgdownload.cse.ucsc.edu/downloads.html
- We download the full data set, but it's possible to interrupt the download during xenoMrna (not needed, too big).
- ``md5sum`` is a basic utility that should be present in your system, otherwise check your packages (yum, apt-get, ...)

.. code-block:: bash

    # location of genome data that can be shared among users
    mkdir -p $GENOMEDIR
    cd $GENOMEDIR
    
    rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/$GENOME/bigZips/ .

    # check downloaded data integrity
    md5sum -c md5sum.txt
    cat *.md5 | md5sum -c

Now unpack the genome - this differs for genomes
some are in single .fa, some are split by chromosomes. Some archives are *tarbombs*, so unpack
to ``chromFa`` directory to avoid possible mess:

.. code-block:: bash    

    mkdir chromFa
    tar xvzf chromFa.tar.gz -C chromFa

Create concatenated genome, use Heng Li's :ref:`sort-alt <sortalt>`
to get common ordering of chromosomes:

.. code-block:: bash

    find chromFa -type f | sort-alt -N | xargs cat > $GENOME.fa

Downlaod all needed annotations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Annotation data is best obtained in UCSC table browser
in BED format and then sorted and indexed by :ref:`BEDtools <bedtools>`

For example: http://genome.ucsc.edu/cgi-bin/hgTables?db=taeGut1:

.. code-block:: bash

    # directory where annotations are stored
    ANNOT=annot
    sortBed -i $ANNOT/ensGene.bed > $ANNOT/ensGene.sorted.bed
    bgzip $ANNOT/ensGene.sorted.bed
    tabix -p bed $ANNOT/ensGene.sorted.bed.gz

FIXME: rozepsat 
Or using compressed files:

.. code-block:: bash

    zcat -d $ANNOT/refSeqGenes.bed.gz | sortBed | bgzip > $ANNOT/refSeqGenes.sorted.bed.gz
    zcat -d $ANNOT/ensGenes.bed.gz | sortBed | bgzip > $ANNOT/ensGenes.sorted.bed.gz
    tabix -p bed $ANNOT/ensGenes.sorted.bed.gz
    tabix -p bed $ANNOT/refSeqGenes.sorted.bed.gz
    
Build indexes for all programs used in the pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Some programs need a preprocessed form of the genome, to speed up their operation.

.. code-block:: bash

    # index chromosome positions in the genome file for samtools
    samtools faidx $GENOMEFA

    # build gmap index for zebra finch
    # with some newer versions it is necessary to use -B <path/to/bindir>
    # beware, this requires quite a lot of memory (gigabytes)
    gmap_build -d $GMAP_IDX -D $GMAP_IDX_DIR $GENOMEFA

    # smalt index
    # needed only for speeding up sim4db
    mkdir -p $GENOMEDIR/smalt
    smalt index -s 4 $SMALT_IDX $GENOMEFA
    
    # convert to blat format
    faToTwoBit $GENOMEFA $GENOME2BIT
    