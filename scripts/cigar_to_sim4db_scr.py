#! /usr/bin/env python

"""Script that parses ``CIGAR`` file 
produced by aligning fragments of contigs to a genome (tested with ``smalt`` output)
and outputs a 'script' for limiting the exon model regions of ``sim4db``.

Input

- output of some read mapper in ``CIGAR`` format

 - all the fragments must be reported by the aligner
 - the fragment names have to be in the same order as the master sequences
 - the fragments must be named like readname_number (like ``fasta_fragments.py`` does)
 - all the hits from one read must follow each other
 
Output

- sim4db 'script'

Algorithm

- load chromosome definition file
- parse the hits:
  - extract read name
  - check if hit is in known chromosome (report error otherwise)
  - for each hit create +-50 KB region clipped to chromosome ends
  - when a read name different from the previous one is encountered, merge all the regions
  - output each contiguos region as a script line

Author: Libor Morkovsky, 2012
"""

# This file is a part of Scrimer.
# See LICENSE.txt for details on licensing.
#    Copyright (C) 2012, 2013 Libor Morkovsky

import sys
import re 

# add a region
def addregion(regs, chr, start, end, reverse):
  # find chromosome + direction combination
  dir = "f" if not reverse else "r"
  
  key = chr + dir
  
  # if this is the first hit of the strand, create a list of regions
  if key not in regs:
    regs[key] = []
  
  # otherwise try to find overlap with some of the existing regions
  regs[key].append((start, end))
  
    
# clip to chromosome length and output the regions  
# merging to overlapping regions on the same chromosome/strand
def output_regions(chromosomes, num, regions):
  """
  from documentation of sim4db:
  A. The input script file format
   
  [-f|-r] -e ESTidx -D GENidx GENlo GENhi
  where:
  cDNAidx        internal index of the cDNA in the input cDNA fasta file
                 (0..#cDNAseqs-1)
  GENidx         internal index of the genomic sequence in the input genome
                 file (0..#GENseqs-1)
  -f             use the cDNA sequence as is
  -r             use the reverse complement of the cDNA sequence
  GENlo, GENhi   begin and end coordinates of the target genomic region;
                 coordinates are 0-based

  For best performance, the script should be sorted by the genomic
  sequence index. 
  """  
  
  for key, ranges in regions.iteritems():
    (seq, dir) = (key[:len(key)-1], key[len(key)-1])
    
    # 
    # merge overlapping
    # 
    ranges.sort()
    newranges = []
    idx = 0
    # go through elements
    while idx < len(ranges):
      # current region
      (start, end) = ranges[idx]
      
      # go through all regions starting after this one
      for nidx in xrange(idx + 1, len(ranges)):
        (nstart, nend) = ranges[nidx]
        # if there is an overlap, try to extend current region
        if nstart <= end:
          end = max(end, nend)
          # continue with merging after this element
          idx = nidx
        # we're in sorted list, there will be no more overlaps
        else:
          break
      
      idx += 1
      newranges.append((start, end))

    for start, end in newranges:
      if start < 0: start = 0
      if end > chromosomes[seq][1]: end = chromosomes[seq][1]
      print "-%s -e %d -D %d %d %d" % (dir, num, chromosomes[seq][0], start, end)

def main():
    if len(sys.argv) != 2:
      exit("use: %s samtools_fai_of_the_genome \n reads cigar from stdin, writes sim4db script to stdout" % sys.argv[0])

    # number of bases to extend on both sides of the hit
    margin = 50000

    # read the fasta index
    fai = open(sys.argv[1])
    chromosomes = dict()
    for num, line in enumerate(fai):
      fields = line.split("\t")
      chromosomes[fields[0]] = (num, int(fields[1]))

    fai.close()

    print >> sys.stderr, len(chromosomes), "sequences loaded from faidx"

    """
     CIGAR format from ssaha2 doc

    [0] in smalt 'cigar' output is preceded by column like 'cigar:S:56'

    [1] query_name  Name of the query sequence from the input file. 
    [2] query_start  Start of the alignment in the query. 
    [3] query_end  End of the alignment in the query. 
    [4] query_strand  Query orientation (+/-). 
    [5] subject_name  Name of the subject sequence from the input file. 
    [6] subject_start Start of the alignment in the subject. 
    [7] subject_end  End of the alignment in the subject. 
    [8] subject_strand  Subject orientation (always + for SSAHA2). 
    [9] score  Smith-Waterman score obtained from the cross_match alignment.
    [10:] series of <operation, length> pairs where operation is one of match, insertion or deletion
    """
    old_seqname = ""
    seqnum = 0
    regions = dict()
    # loop through the input
    for nline, line in enumerate(sys.stdin):
      fields = line.split()
      seqname = re.sub("_[0-9]+$", "", fields[1])
      
      # if new query sequence encountered, output the regions
      if seqname != old_seqname and nline != 0:
        if len(regions):
            output_regions(chromosomes, seqnum, regions)
            regions = dict()

        seqnum += 1
      
      # store the seqname for next loop
      old_seqname = seqname
      
      # sanity check -- if we're matching against the correct genome
      if fields[5] not in chromosomes:
        # a special case for 'no hit'
        if fields[5] != '*':
            print >> sys.stderr, "hit %d was matched to sequence not in faidx/genome (%s)" % (nline, fields[5])
        continue
      
      # create a region for current hit
      # subj_name, start, end, reverse?
      addregion(regions, fields[5], int(fields[6]) - margin, int(fields[7]) + margin, fields[4] != fields[8])

    # output the ranges for last sequence
    output_regions(chromosomes, seqnum, regions)

if __name__ == '__main__': main()
