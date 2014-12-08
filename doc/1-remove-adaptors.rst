Remove cDNA synthesis adaptors
==============================

Quality check of the raw data
-----------------------------
Results of the quality check of the raw data can be used as a reference point
to check the improvements done by this step.

.. code-block:: bash

    IN=00-raw
    OUT=10-fastqc
    mkdir $OUT
    fastqc --outdir=$OUT --noextract --threads=$CPUS $IN/*.fastq

Split the files by multiplex tags
---------------------------------
If you used multiplexing during library preparation, your reads have either
already been split in your sequencing facility, or you have to split the 
reads according to the 'barcodes' stored in the sequences yourself. 

- for 454, this is done by ``sfffile``, you have to convert resulting sff files to fastq format
- for Illumina this can be done e.g. by `eautils <https://code.google.com/p/ea-utils/>`_

.. code-block:: bash

    IN=03-sff
    OUT=04-sff-split
    mkdir -p $OUT

    parallel -j 1 sfffile -s RLMIDs -o $OUT/{/.} {} ::: $IN/*.sff

Remove cDNA synthesis primers with cutadapt
-------------------------------------------
Another source of noise in the data is the primers that were used for reverse transcription
of mRNA and for the following PCR amplification of cDNA. We remove them using ``cutadapt``.

If you have more or less than three primers, you have to change the command by adding/removing
``--anywhere=`` parts.

.. code-block:: bash
    
    # data from previous steps
    IN=10-mid-split
    OUT=12-cutadapt
    mkdir $OUT

    # cut out the evrogen sequences using GNU parallel and cutadapt
    # cutadapt supports only 'N' wildcards, no ambiguity codes
    parallel -j $CPUS "cutadapt --anywhere=$PRIMER1 --anywhere=$PRIMER2 --anywhere=$PRIMER3 \
      --error-rate=0.2 --overlap=15 --minimum-length=40 \
      --output=$OUT/{/.}.fastq --rest-file=$OUT/{/.}.rest {}" ::: $IN/*.fastq > $OUT/cutadapt.log

Check the results
-----------------
It is necessary to check the results of adaptor cutting. 

First you can check how many of the primers were missed by cutadapt. ``agrep`` uses a different 
matching algorithm than cutadapt, so some remaining hits are usually found.
``/dev/null`` is used as the second input to ``agrep`` so the filenames are output.

.. code-block:: bash

    NERR=5
    parallel agrep -c -$NERR "$PRIMER3" {} /dev/null ::: $OUT/*.fastq|grep -v /dev

The next thing to check is the logs produced by ``cutadapt``.
Results for our data - 454 Titanium data from Smart kit synthesized cDNA: 

- ~70% trimmed reads
- ~10% trimmed basepairs
- ~10% too short reads

.. code-block:: bash

    grep -A5 Processed $OUT/cutadapt.log | less

The mean length of the removed sequences should be close to length of the adapter (31 in this case):

.. code-block:: bash

    less $OUT/cutadapt.log

.. code-block:: bash

    # Lengths of removed sequences (5')
    # length  count   expected
    # 5       350     264.7
    # 6       146     66.2
    # ...
    # 30      6414    0.0
    # 31      63398   0.0
    # 32      6656    0.0
    # ...

The size of the ``.rest`` files is 1/500 of the ``.fastq`` (should be 1/250 for ``.fasta``)

.. code-block:: bash

    ls -l $OUT

The ``fastqc`` checks should be +- ok.

.. code-block:: bash

    fastqc --outdir=13-fastqc --noextract --threads=8 $OUT/*.fastq

Visual debugging
----------------
If something in the previous checks looks weird, look directly at the data. Substitute filenames below with 
the names of your files. 

Look where the primers are in the sequence. ``tre-agrep`` is used to color the output of ``agrep``, because
``agrep`` throughput is ~ 42 MB/s while ``tre-agrep`` throughput is ~ 2 MB/s.

.. code-block:: bash

    FQFILE=$IN/G3UKN3Q01.fasta
    NERR=5
    agrep -n -$NERR "$PRIMER3" $FQFILE |tre-agrep -$NERR "$PRIMER3" --color|less -S -R

To find out how many differences should be allowed in the pattern matching, we try to find a value of ``NERR``
where the primer sequence starts to match randomly inside the reads, and not only in at the beginning.
Notice the ``^`` in the first command, marking the start of the read.

.. code-block:: bash

    agrep -c -$NERR "^$PRIMER3" $FQFILE && agrep -c -$NERR "$PRIMER3" $FQFILE

    # numbers for tag-cleaned G59B..
    # 4 errors: 11971 12767
    # 5 errors: 16366 17566
    # 6 errors: 17146 23858
    # 7 errors: 18041 67844

In our sample results, numbers start to diverge for ``NERR`` > 5, so 5 is a good choice.

Read count statistics
---------------------

For a single file:

.. code-block:: bash

    # read count statistics
    # @ can be in the beginning of quality string, so filter the rows in order

    # count of sequences
    awk '((NR%4)  == 1)' $FQFILE | wc -l
    # or more effective
    echo $(( $(wc -l $FQFILE) / 4 ))

    # count of sequenced bases
    awk '((NR%4)  == 2)' $FQFILE | wc -m

For all files in ``OUT``:

.. code-block:: bash

    # parallel, IO bound task, so run one process a time
    OUT=12-cutadapt
    echo "read_count base_count filename"
    parallel -j 1 'echo $( gawk "((NR%4)  == 1)" {} | wc -l ) $( gawk "((NR%4)  == 2)" {} | wc -m ) {}' ::: $OUT/*.fastq
