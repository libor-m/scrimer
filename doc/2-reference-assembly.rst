Assemble reads into contigs
===========================

Use newbler to assemble the reads
---------------------------------
Here **newbler** is used to assemble the contigs. For 3 GB of read data the assembly took 25 CPU hours and 15 GB RAM.

.. code-block:: bash

    IN=12-cutadapt
    OUT=22-newbler
    CPUS=19

    # run full assembly 
    runAssembly -o $OUT -cdna -cpu $CPUS -m $IN/*.fastq


Remove contigs that are similar to each other
---------------------------------------------
The aim is to get one transcript per locus, preferably the longest one. Otherwise the read mapping
process would be faced with many ambiguous locations. We achieve this by:

- doing all-to-all comparison within the isotigs
- grouping isotigs that are similar up to given thresholds of coverage and identity,
  (disjoint sets forest graph algorithm is used)
- and selecting only the longest contig for each group

.. code-block:: bash

    IN=22-newbler
    OUT=$IN
    SEQFILE=$IN/454Isotigs.fna 
    MINCOVERAGE=90
    MINIDENTITY=95

    # taken from lastz human-chimp example, should be report only highly similar hits
    # filter self matches with awk
    lastz $SEQFILE[multiple] $SEQFILE \
          --step=10 --seed=match12 --notransition --exact=20 --noytrim \
          --match=1,5 --ambiguous=n \
          --coverage=$MINCOVERAGE --identity=$MINIDENTITY \
          --format=general:name1,size1,start1,name2,size2,start2,strand2,identity,coverage \
          | awk '($1 != $4)' > $SEQFILE.lastz-self

    # takes 8 minutes, finds 190K pairs

    # find the redundant sequences
    tail -n +2 $SEQFILE.lastz-self | find_redundant_sequences.py > $SEQFILE.redundant

    # add the short sequences to discard list
    grep '>' $SEQFILE | awk '{ sub(/length=/,"",$3); sub(/^>/, "", $1); if($3 - 0 < 300) print $1;}' >> $SEQFILE.redundant

    # get rid of the redundant ones
    seq_filter_by_id.py $SEQFILE.redundant 1 $SEQFILE fasta - $SEQFILE.filtered


Check the results
-----------------

.. code-block:: bash

    # check the input contig length distribution
    grep '>' $SEQFILE | awk '{ sub(/length=/,"",$3); print $3}' | histogram.py -b 30

    # view how many contigs we got after filtering
    grep -c '>' $SEQFILE.filtered

Assembling Illumina data
------------------------

`Trinity <http://trinityrnaseq.sourceforge.net/>`_ gives fairly good results for transcriptome assembly from Illumina data.
Simple Trinity usage looks like:

.. code-block:: bash

    # as mentioned on the homepage
    Trinity --seqType fq --JM 50G --left reads_1.fq  --right reads_2.fq --CPU 6

Some more `tips on assembling 'perfect' transcripts <ftp://flamingo.bio.indiana.edu/evigene/docs/perfect-mrna-assembly-2013jan.txt>`_ by Don Gilbert.
