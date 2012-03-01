#!/usr/bin/env python

#
# per base statistics of 5' ends of given fastq input
#  useful for finding/recognizing primer sequences in the reads
#
# use: 5prime-stats.py fastq [n-bases=30]

from Bio import SeqIO
import sys

file = SeqIO.parse(sys.argv[1], "fastq")

# stats for first n chars
nbases = 30
if len(sys.argv) > 1:
  nbases = int(sys.argv[2])
d = [dict() for x in range(nbases)]

for rec in file:
  for i in range(min(len(d), len(str(rec.seq)))):
    letter = str(rec.seq)[i]
    if letter not in d[i]:
      d[i][letter] = 0
      
    d[i][letter] += 1

# output the stats
from numpy import var, array
topseq = ""
for i in range(len(d)):
  # normalize the histogram to the frequencies
  res = array([float(x)/sum(d[i].values()) for x in d[i].values()])
  # print res
  idx = res.argmax()
  topletter = d[i].keys()[idx]
  print i, ": ", topletter, res[idx], var(res)
  topseq += topletter
  
print "top sequence:", topseq


