#!/usr/bin/env python

# given a fastq file and a list of sequence indices to remove
# outputs the file without the given sequences (index starts at 0)
#
# this could probably be simple sed, if it didn't hit the command line size limit

import sys

# read the kill list
fkill = open(sys.argv[1])
killlines = fkill.readlines()
killseq = frozenset([int(line) for line in killlines])
fkill.close()

# go through the input, avoid outputting the 4 line blocks 
# with given indices
fseq = open(sys.argv[2])
linenum = 0
for line in fseq:
    seqnum = linenum / 4
    linenum += 1  
    if seqnum in killseq:
        continue
    print line.strip()

fseq.close()