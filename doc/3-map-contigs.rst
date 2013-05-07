Map contigs to the reference genome
===================================

Mapping options
---------------

GMAP
^^^^
.. code-block:: bash


    INFILE=20-jp-contigs/lu_master500_v2.fna.filtered
    OUT=30-tg-gmap
    mkdir $OUT
    
    OUTFILE=$OUT/lu_master300_v2.gmap.gff3

    gmap -D $GMAP_IDX_DIR -d $GMAP_IDX -f gff3_gene -B 3 -x 30 -t $CPUS\
        --cross-species $INFILE  > $OUTFILE


sim4db
^^^^^^
Sim4db is considerably slow, but it has an option to provide a *mapping script*. We exploit this by 
creating a mapping script by using fast short read aligner with fragments of contigs we want to map,
and then convert hits of those short reads into candidate mapping locations.

Please use the patched version sim4db to obtain correct mapping coordinates in the gff files.

.. code-block:: bash

    INFILE=20-jp-contigs/lu_master500_v2.fna.filtered
    OUT=31-tg-sim4db
    mkdir $OUT

    # these values are derived, it's not necessary to change them
    FRAGS=$OUT/${INFILE##*/}.frags
    SMALT_OUT=$FRAGS.cigar
    SIM4_SCR=${FRAGS%.*}.sim4scr
    OUT0=${FRAGS%.*}.tmp.gff3
    OUTFILE=${FRAGS%.*}.gff3

Use ``smalt`` as a fast mapper to find all +-50 kBase windows for predicting 
exon/gene models with sim4db:

.. code-block:: bash

    # create fragments, using slightly modified fasta_fragments.py from lastz distribution
    cat $INFILE | fasta_fragments.py --step=80 > $FRAGS

    # map the fragments with smalt (takes few minutes), reporting all hits (-d -1) scoring over 60
    smalt_x86_64 map -n 8 -f cigar -o $SMALT_OUT -d -1 -m 60 $SMALT_IDX $FRAGS

    # construct the script for sim4db
    cat $SMALT_OUT | cigar_to_sim4db_scr.py $GENOMEFA.fai | sort --key=5n,5 > $SIM4_SCR

Run sim4db using the script. (takes several seconds for the whole genome):

.. code-block:: bash

    sim4db -genomic $GENOMEFA -cdna $INFILE -script $SIM4_SCR -output $OUT0 -gff3 -interspecies -mincoverage 70 -minidentity 90 -minlength 60 -alignments -threads $CPUS

    # fix chromosome names 
    sed s/^[0-9][0-9]*:chr/chr/ $OUT0 > $OUTFILE

Transfer genome annotations to our contigs
------------------------------------------
Annotatate our sequences by data from similar sequences in reference genome.

    sim4db manual [#]_ states: Exon coordinates are nucleotide based, starting from 1. Genomic coordinates are always 
    in the original sequence, while the cDNA coordinates will refer to positions in the reverse 
    complement of the sequence if the match orientation is indicated as 'complement'.

- this seems unnecessary, because the orientation of the transcript can be deduced from the target chromosome strand (..?)
- patched in sim4db, submitted patch to sourceforge

Each contig mapping to genome creates different coordinate system for transferring
the annotations. Annotation data should be sorted and tabix indexed. Multiple coordinate systems 
and multiple annotations can be used. There is one output per coordinate system, they're merged 
during the scaffold phase:

.. code-block:: bash

    # multiple coordinate systems if needed (one system per mapping)
    COORDS="30-tg-gmap/lu_master300_v2.gmap.gff3 31-tg-sim4db/lu_master500_v2.fna.filtered.gff3"
    ANNOTS=/data/genomes/taeGut1/annot/ensGene_s.bed.gz
    OUT=32-liftover
    mkdir -p $OUT

    for C in $COORDS
    do
        liftover.py "$C" $ANNOTS > $OUT/${C##*/}-lo.gff3
    done  


Create 'transcript scaffold' using the annotations
--------------------------------------------------
Construct a 'transcript scaffold' (contigs joined in order of appearance on reference genome chromosomes).
This is mainly because of viewing conveninence with IGV. 'N' gaps should be larger than max read size
to avoid the mapping of the reads across gaps:

.. code-block:: bash
    
    # filtered contigs
    INFILE=20-jp-contigs/lu_master500_v2.fna.filtered
    # transferred annotations from previous step
    ANNOTS=32-liftover/*-lo.gff3
    # output directory
    OUT=33-scaffold
    # name of the output 'genome'
    GNAME=lx4

    mkdir $OUT
    OUTGFF=$OUT/$GNAME.gff3

    scaffold.py $INFILE $ANNOTS $OUT/$GNAME.fasta $OUTGFF

    # sort, compress and index the merged annotations
    # so they can be used further down in the pipeline
    OUTFILE=${OUTGFF%.*}.sorted.gff3

    sortBed -i $OUTGFF > $OUTFILE
    bgzip $OUTFILE
    tabix -p gff $OUTFILE.gz

.. [#] http://sourceforge.net/apps/mediawiki/kmer/index.php?title=Getting_Started_with_Sim4db
