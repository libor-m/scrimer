#------------------------------------------------
# Operation
#------------------------------------------------

#  1. split the files according to MIDs with SFFile
#  2. remove cDNA synthesis primers with cutadapt
#  3. check the results

#  1. split the files according to MIDs with SFFile
#------------------------------------------------

#
# sff format is necessary for this and we don't have any
#

#  2. remove cDNA synthesis primers with cutadapt
#------------------------------------------------

###
# primers used for library construction (taken from evrogen .pdf report)
# AAGCAGTGGTATCAACGCAGAGTACGCGGG
# AAGCAGTGGTATCAACGCAGAGT
# AAGCAGTGGTATCAACGCAGAGTTTTTGTTTTTTTCTTTTTTTTTTVN
###

IN=10-mid-split
OUT=12-cutadapt
NERR=5
mkdir $OUT

# cut out the evrogen sequences using GNU parallel and cutadapt
# cutadapt supports only 'N' wildcards, no ambiguity codes
parallel cutadapt --anywhere=AAGCAGTGGTATCAACGCAGAGTTTTTGTTTTTTTCTTTTTTTTTTNN --anywhere=AAGCAGTGGTATCAACGCAGAGTACGCGGG --anywhere=AAGCAGTGGTATCAACGCAGAGT \
--error-rate=0.2 --overlap=15 --minimum-length=40 \
--output=$OUT/{/.}.fastq --rest-file=$OUT/{/.}.rest {} ::: $IN/*.fastq > $OUT/cutadapt.log

#  3. check the results
#------------------------------------------------

# check the number of remaining hits (using the /dev/null trick to get the filenames and filtering it out by grep)
parallel agrep -c -$NERR "AAGCAGTGGTATCAACGCAGAGT" {} /dev/null ::: $OUT/*.fastq|grep -v /dev

# look at the log, check the sanity of output
grep -A5 Processed $OUT/cutadapt.log | less
# results for 454 Titanium data from Smart kit synthesized cDNA, luscinia: 
#  ~70% trimmed reads
#  ~10% trimmed basepairs
#  ~10% too short reads

less $OUT/cutadapt.log
# length of the removed sequence should be equal to length of the adapter (31 in this case):

# Lengths of removed sequences (5')
# length  count   expected
# 5       350     264.7
# 6       146     66.2
# ...
# 30      6414    0.0
# 31      63398   0.0
# 32      6656    0.0
# ...

ls -l $OUT
# size of the .rest files is 1/500 of the .fastq (should be 1/250 for .fasta)

# now the fastqc checks should be +- ok
fastqc --outdir=13-fastqc --noextract --threads=8 $OUT/*.fastq

#------------------------------------------------
# Visualization
#------------------------------------------------

# if something went wrong, look at the data directly
FQFILE=$IN/G3UKN3Q01.fasta

# check the results, use the fast agrep for searching, (really) slow tre-agrep 
# for coloring of the approximate matches
# use -n to assess how often does a tag occur
agrep -n -$NERR "AAGCAGTGGTATCAACGCAGAGT" $FQFILE |tre-agrep -$NERR "AAGCAGTGGTATCAACGCAGAGT" --color|less -S -R

# finding an 'elbow' in error tolerace - difference between count of ^TAG and TAG matches
# until the numbers are close, the number of allowed errors is ok
# this exploits the fact that most of the real tags is in the beginning of the sequences
# beware, this does not work for systematic error in the files, like other 4 bases prepended 
# in to the seqnece for NERR=5 (it fits into the tolerance)
agrep -c -$NERR "^AAGCAGTGGTATCAACGCAGAGT" $FQFILE && agrep -c -$NERR "AAGCAGTGGTATCAACGCAGAGT" $FQFILE
# numbers for tag-cleaned G59B..
# 4 errors: 11971 12767
# 5 errors: 16366 17566
# 6 errors: 17146 23858
# 7 errors: 18041 67844
