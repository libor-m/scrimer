Prepare reference genome
========================

Downlaod and prepare the reference genome
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- List of available genomes is at http://hgdownload.cse.ucsc.edu/downloads.html
- We download the full data set, but it's possible to interrupt the download during xenoMrna (not needed, too big).
- ``md5sum`` is a basic utility that should be present in your system, otherwise check your packages (yum, apt-get, ...)

:: 

    # location of genome data that can be shared among users
    cd $GENOMES
    mkdir $GENOME
    cd $GENOME
    rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/$GENOME/bigZips/ .

    # check downloaded data integrity
    md5sum -c md5sum.txt
    cat *.md5|md5sum -c

Unpack the genome - this differs for genomes
some are in single .fa, some are split by chromosomes ::
    
    tar xvzf chromFa.tar.gz

Create concatenated genome, use Heng Li's sort-alt
to get common ordering of chromosomes::

    find chromFa -type f|sort-alt -N|xargs cat > $GENOME.fa

Downlaod all needed annotations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Annotation data is best obtained in UCSC table browser
in BED format and then sorted and indexed by BEDtools

For example: http://genome.ucsc.edu/cgi-bin/hgTables?db=taeGut1::

    # directory where annotations are stored
    ANNOT=annot
    sortBed -i $ANNOT/ensGene.bed > $ANNOT/ensGene.sorted.bed
    bgzip $ANNOT/ensGene.sorted.bed
    tabix -p bed $ANNOT/ensGene.sorted.bed.gz

Build indexes for all programs used in the pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Some programs need a preprocessed form of the genome, to speed up their operation.

::

    # index chromosome positions in the genome file for samtools
    samtools faidx $GENOMES/$GENOME/$GENOME.fa

    # build gmap index for zebra finch
    gmap_build -d gmap_taeGut1 -D $GENOMES/$GENOME $GENOMES/$GENOME/$GENOME.fa


    # smalt index
    # recommended settings for 454 (step 4, k-mer size 13)
    mkdir smalt
    smalt index -s 4 smalt/${GENOME}k13s4 $GENOME.fa
    