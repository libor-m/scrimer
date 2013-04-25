#!/usr/bin/env python

"""
Find the most common letter in first n bases of reads in FASTQ file.
Useful for finding and recognizing primer sequences in the reads.
"""

# This file is a part of Scrimer.
# See LICENSE.txt for details on licensing.
#    Copyright (C) 2012, 2013 Libor Morkovsky

from Bio import SeqIO
import sys

def main():
    if len(sys.argv) < 2:
        print >> sys.stderr, "use: %s fastq_file [n-bases=30]" % sys.argv[0]

    file = SeqIO.parse(sys.argv[1], "fastq")

    # stats for first n chars
    nbases = 30
    if len(sys.argv) > 2:
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

if __name__ == '__main__': main()