#
# http://biostar.stackexchange.com/questions/3468/how-to-best-deal-with-adapter-contamination-illumina
#
# find out about the primers (see file clontech_primers)
export FQFILE=G59B7NP01.fastq
# count of all sequences
grep -c "^@G59B" G59B7NP01.fastq
# 756882

# sequences having 5' pcr primer - or cDNA synth primers in the beginning
grep -c "^AAGCAGTGGTATCAACGCAGAGT" $FQFILE
# 387224

# sequences containing the primer sequence
grep -c "AAGCAGTGGTATCAACGCAGAGT" $FQFILE
# 390224

echo $(( 390224 - 387224 ))
# 3000 sequences not having the sequence at beginning, let's look at it
grep -v "^AAGCAGTGGTATCAACGCAGAGT" $FQFILE | grep "AAGCAGTGGTATCAACGCAGAGT" --color=always|less -S -R
# many of those sequences look a bit chimeric, the rest has one more A at the beginning

# count the chimeric ones
grep -v '^.\?AAGCAGTGGTATCAACGCAGAGT' $FQFILE | grep -c "AAGCAGTGGTATCAACGCAGAGT"
grep -v '^.\?AAGCAGTGGTATCAACGCAGAGT' $FQFILE | grep  "AAGCAGTGGTATCAACGCAGAGT" --color=always|less -S -R
# 1519

# count the SMARTer V oligos
agrep -1 -c "AAGCAGTGGTATCAACGCAGAGTACGCGGGG" $FQFILE
# 188159
# seems like it's worth removing those first, and then the PCR primer

# good option is to discard such sequences
# count the discarded sequences with approximate matching
agrep -v -2 "^AAGCAGTGGTATCAACGCAGAGT" $FQFILE | agrep -2 -c "AAGCAGTGGTATCAACGCAGAGT"
# 2060 sequences, not a big deal..

# any reverse complements?
agrep -1 -c $SEQ $FQFILE
SEQ=
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
ctech_5prime = Seq("AAGCAGTGGTATCAACGCAGAGT", IUPAC.unambiguous_dna)
ctech_3prime = Seq("CGGGGTACGATGAGACACCA", IUPAC.unambiguous_dna)
ctech_5prime.reverse_complement()
# Seq('ACTCTGCGTTGATACCACTGCTT', IUPACUnambiguousDNA()) -> 1322 hits
ctech_5prime.complement()
# Seq('TTCGTCACCATAGTTGCGTCTCA', IUPACUnambiguousDNA()) -> 0 hits
 "AAGCAGTGGTATCAACGCAGAGT"[::-1]
# 'TGAGACGCAACTATGGTGACGAA' -> 0 hits
ctech_3prime.reverse_complement()
# Seq('TGGTGTCTCATCGTACCCCG', IUPACUnambiguousDNA()) -> 0 hits
ctech_3prime.complement()
# Seq('GCCCCATGCTACTCTGTGGT', IUPACUnambiguousDNA()) -> 0 hits
'CGGGGTACGATGAGACACCA'[::-1]
# 'ACCACAGAGTAGCATGGGGC' (reverse of ctech_3prime) -> 0 hits

# why there are no hits form the other PCR primer?
# why there are so many hits (50%) with the 5' primer in the beginning?
#  i'd expect more of in-the-middle fragments present in the nebulized cDNA...?

#
# let's remove the tags..
# firts the SMARTer V, then the PCR primer
tagcleaner.pl -fastq 10-split/G59B7NP01.fastq -out 12-remove-primers/G59B7NP01-a -minlen 50 -tag5 AAGCAGTGGTATCAACGCAGAGTACGCGGG -mm5 2 -split
tagcleaner.pl -fastq 12-remove-primers/G59B7NP01-a.fastq -out 12-remove-primers/G59B7NP01 -minlen 50 -tag5 AAGCAGTGGTATCAACGCAGAGT -mm5 2 -split


tagcleaner.pl -fastq 10-split/GSVZDOM02.fastq -out 12-remove-primers/GSVZDOM02-a -minlen 50 -tag5 AAGCAGTGGTATCAACGCAGAGTACGCGGG -mm5 2 -split
tagcleaner.pl -fastq 12-remove-primers/GSVZDOM02-a.fastq -out 12-remove-primers/GSVZDOM02 -minlen 50 -tag5 AAGCAGTGGTATCAACGCAGAGT -mm5 2 -split

# check the results, use the fast agrep for searching, (really) slow tre-agrep 
# for coloring of the approximate matches
export NERR=5 
agrep -$NERR "AAGCAGTGGTATCAACGCAGAGT" $FQFILE |tre-agrep -$NERR "AAGCAGTGGTATCAACGCAGAGT" --color|less -S -R
# 621 hits in the G5.. file

# finding an 'elbow' in error tolerace - difference between count of ^TAG and TAG matches
# until the numbers are close, the number of allowed errors is ok
# this exploits the fact that most of the real tags is in the beginning of the sequences
agrep -c -$NERR "^AAGCAGTGGTATCAACGCAGAGT" $FQFILE && agrep -c -$NERR "AAGCAGTGGTATCAACGCAGAGT" $FQFILE
# numbers for tag-cleaned G59B..
# 4 errors: 11971 12767
# 5 errors: 16366 17566
# 6 errors: 17146 23858
# 7 errors: 18041 67844

# there's no easy way to use agrep directly for the tag filtering/splitting
# - it does not output the real sequence matched (-o in grep) nor offsets of the match ;(
# the conservative resolution for now is to remove all the matching sequences
# construct the killlist (each record in fastq has 4 lines)

# beware! agrep messes up the line numbers if lines are longer than 1024 chars
agrep -un -5 "AAGCAGTGGTATCAACGCAGAGT" $FQFILE|tr -d ':'|gawk '{print int($1/4);}' > $FQFILE.kill

# use float arithmetic in gawk as a check (then issue grep -c '\.' $FQFILE.kill)
# (sequence is always in second line, we want 0 based index as output, 4 should divide each line num)
# pv is just for monitoring of the slow tre- search..;)
tre-agrep -n -5 "AAGCAGTGGTATCAACGCAGAGT" $FQFILE|cut -d: -f1|gawk '{print ($1-2)/4;}'|pv > $FQFILE.kill
./fastq_kill_lines.py $FQFILE.kill $FQFILE > $FQFILE.new
mv $FQFILE.new $FQFILE && rm $FQFILE.kill

# TODO: rescue the ^TAG reads with misread tags

#
# another tag cleaner tool: cutadapt, uses gapped alignments, in python
#
# given 2 sequences, only the best matching should be trimmed
# error rate of 0.2 should be equivalent to 5 mismatches out of 23 base pair tag
cutadapt --anywhere=AAGCAGTGGTATCAACGCAGAGTACGCGGGG --anywhere=AAGCAGTGGTATCAACGCAGAGT --error-rate=0.2 --overlap=5 --minimum-length=40 --output=13-cutadapt/${FQFILE##*/} --rest-file=13-cutadapt/${FQFILE##*/}.rest $FQFILE
# in a loop 
for FQFILE in 10-split/*.fastq; do cutadapt --anywhere=AAGCAGTGGTATCAACGCAGAGTACGCGGGG --anywhere=AAGCAGTGGTATCAACGCAGAGT --error-rate=0.2 --overlap=5 --minimum-length=40 --output=13-cutadapt/${FQFILE##*/} --rest-file=13-cutadapt/${FQFILE##*/}.rest $FQFILE ; done 2>&1 > 13-cutadapt/cutadapt.log

agrep -c -$NERR "AAGCAGTGGTATCAACGCAGAGT"  13-cutadapt/*.fastq
# seems like only a bearable number of sequences contain the pcr primer, ignore it for now
13-cutadapt/G59B7NP01.fastq: 0 # this one was fastq_kill.. cleaned
13-cutadapt/G60Z2EH01.RL11.fastq: 224
13-cutadapt/G60Z2EH01.RL9.fastq: 310
13-cutadapt/G60Z2EH02.RL10.fastq: 369
13-cutadapt/G60Z2EH02.RL12.fastq: 268
13-cutadapt/GS60IET02.RL1.fastq: 299
13-cutadapt/GS60IET02.RL2.fastq: 286
13-cutadapt/GSVZDOM01.RL1.fastq: 8
13-cutadapt/GSVZDOM01.RL2.fastq: 25
13-cutadapt/GSVZDOM02.fastq: 96
13-cutadapt/GYAB93P02.RL1.fastq: 327
13-cutadapt/GYAB93P02.RL2.fastq: 408

###
# primers really used for library construction (taken from .pdf report)
# AAGCAGTGGTATCAACGCAGAGTACGCGGG
# AAGCAGTGGTATCAACGCAGAGT
# AAGCAGTGGTATCAACGCAGAGTTTTTGTTTTTTTCTTTTTTTTTTVN
###

# cut out the evrogen sequences using GNU parallel and cutadapt
# cutadapt supports only 'N' wildcards, no ambiguity codes
parallel cutadapt --anywhere=AAGCAGTGGTATCAACGCAGAGTTTTTGTTTTTTTCTTTTTTTTTTNN --anywhere=AAGCAGTGGTATCAACGCAGAGTACGCGGG --anywhere=AAGCAGTGGTATCAACGCAGAGT \
--error-rate=0.2 --overlap=5 --minimum-length=40 \
--output=14-cutadapt-2/{/.}.fastq --rest-file=14-cutadapt-2/{/.}.rest {} ::: 10-split/*.fastq > 14-cutadapt-2/cutadapt.log

# check the number of remaining hits (using the /dev/null trick to get the filenames and filtering it out by grep)
parallel agrep -c -5 "AAGCAGTGGTATCAACGCAGAGT" {} /dev/null ::: 14-cutadapt-2/*.fastq|grep -v /dev
14-cutadapt-2/GSVZDOM01.RL1.fastq: 7
14-cutadapt-2/GSVZDOM01.RL2.fastq: 23
14-cutadapt-2/GSVZDOM02.fastq: 95
14-cutadapt-2/G60Z2EH01.RL11.fastq: 217
14-cutadapt-2/G60Z2EH02.RL12.fastq: 264
14-cutadapt-2/G60Z2EH01.RL9.fastq: 299
14-cutadapt-2/GS60IET02.RL1.fastq: 285
14-cutadapt-2/G60Z2EH02.RL10.fastq: 341
14-cutadapt-2/GS60IET02.RL2.fastq: 274
14-cutadapt-2/GYAB93P02.RL1.fastq: 313
14-cutadapt-2/GYAB93P02.RL2.fastq: 386
14-cutadapt-2/G59B7NP01.fastq: 814

