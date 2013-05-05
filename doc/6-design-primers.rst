Design primers
==============
.. note::
    
    Ideal case would be to add the annotations found by ``isPcr`` to the primer gff3 tags.

Design primers with Primer3
---------------------------
Set inputs and outputs for this step::

    VARIANTS=50-variants/lx4-variants-selected.vcf.gz
    CONTIGS=33-scaffold/lx4.fasta
    ANNOTS=33-scaffold/lx4.sorted.gff3.gz 

    # other inputs/outputs
    GFFILE=lx4-primers.gff3
    OUT=60-gff-primers
    GFF=$OUT/$GFFILE
    mkdir -p $OUT

    # for all selected variants design pcr and genotyping primers
    # takes about a minute for 1000 selected variants, 5 MB gzipped vcf, 26 MB uncompressed genome, 5 MB gzipped gff
    design_primers.py $CONTIGS $ANNOTS $VARIANTS > $GFF

It's quite good to sort and index the annotation before using it in IGV::

    sortBed -i $GFF | bgzip > ${GFF%.*}.sorted.gff3.gz
    tabix -f -p gff ${GFF%.*}.sorted.gff3.gz

Validate primers with isPcr
---------------------------
# first, check primers with sole blat
# the next step is to set up a gfServer and use gfPcr
# extract the primers
./extract_primers.py 60-gff-primers/lx3-primers2.sorted.gff3.gz isPcr > 61-check-primers/lx3-primers2.isPcr
./extract_primers.py 60-gff-primers/lx3-primers2.sorted.gff3.gz > 61-check-primers/lx3-primers2.fa

# check against transcriptome data
isPcr -out=psl 61-check-primers/lx3.2bit 61-check-primers/lx3-primers2.isPcr 61-check-primers/lx3-primers2.isPcr.lx3.psl
blat -minScore=15 -tileSize=11 -stepSize=5 -maxIntron=0 61-check-primers/lx3.2bit 61-check-primers/lx3-primers2.fa 61-check-primers/lx3-primers2.lx3.psl

# check against reference genome
isPcr -out=psl /data/genomes/taeGut1/twoBit/taeGut1.2bit 61-check-primers/lx3-primers2.isPcr 61-check-primers/lx3-primers2.isPcr.taeGut1.psl
blat -minScore=15 -tileSize=11 -stepSize=5 -maxIntron=0 /data/genomes/taeGut1/twoBit/taeGut1.2bit 61-check-primers/lx3-primers2.fa 61-check-primers/lx3-primers2.taeGut1.psl

###
# use blat and isPcr to map previously manually designed primers
###

# convert fa to 2bit
faToTwoBit 33-virtual-genome/lx3.fasta 61-check-primers/lx3.2bit

# use isPcr to check the products
isPcr 61-check-primers/lx3.2bit 61-check-primers/manual_pcr.ispcr 61-check-primers/manual_pcr.ispcr.fa

# isPcr with psl output to be easily loaded into IGV
isPcr -out=psl 61-check-primers/lx3.2bit 61-check-primers/manual_pcr.ispcr 61-check-primers/manual_pcr.ispcr.psl

# blat the manually designed primers againts the virtual genome
# primer sequences are short, so the blat default parameters
# have to be changed a bit
blat -minScore=15 -tileSize=6 -maxIntron=0 61-check-primers/lx3.2bit 61-check-primers/manual_pcr.fa 61-check-primers/manual_pcr.psl

# pcr recomended parameters for blat:  tileSize=11, stepSize=5
# http://genomewiki.ucsc.edu/index.php/Blat-FAQ
blat -minScore=15 -tileSize=11 -stepSize=5 -maxIntron=0 61-check-primers/lx3.2bit 61-check-primers/manual_pcr.fa 61-check-primers/manual_pcr-2.psl

# agrep is quite enough for simple checks on assemblies of this size (30 MB)
SEQ=GCACATTTCATGGTCTCCAA
agrep $SEQ 0a-jp-newbler-contigs/lu??_contigs.fasta|grep $SEQ


Visualization
=============
See all places where ``primer3`` reported problems::

    grep primer-gt $GFF | grep -v -c 'PROBLEMS='
