#! /usr/bin/env python

# input
#  fasta file with the original sequences
#  set of gff files with exon features having the Target attribute 
#   (product of liftover.py)
#  other gff/bed files to remap to the genome
#
# output 
#  fasta with the input sequences pasted with 50 Ns
#   - unambiguos sequences assigned to 'chromosomes' in the order 
#   of the template genome (that was used to generate the Target exon mappings)
#   - ambiguos (having more than one possible chromosome) assigned to chrAmb
#   - unmapped sequences assigned to chrUnVirt (to distinguish from chrUn in target genome)
#  gff file with the locations of the input sequences (gene)
#   and remapped contents of the input gff files
# 
# algorithm
#  go through the input gff files, construct a dictionary {read_name -> {chr -> location}}
#  add the lowest coordinate found on a given chromosome (overwriting previous values)
#  sort the list with single candidate locations with 'chromocompare'
# 
# Author: Libor Morkovsky 2012
#

import sys
import pybedtools

# reuse the generator fasta reader from lastz
from fasta_fragments import fasta_sequences

def main():
  if len(sys.argv) < 2:
    sys.exit('use: %s input_sequences gff_with_mapped_exons [...]' % sys.argv[0])
  
  # load the sequences into dict {read_name -> sequence}
  seqfile = open(sys.argv[1])
  reads = { read_name : seq for (read_name, seq) in fasta_sequences(seqfile) }
  seqfile.close()
  
  # construct {read_name -> {chr -> location}}
  readloc = dict()
  features = dict()
  gffiles = [pybedtools.BedTool(fname) for fname in sys.argv[2:]]
  for gff in gffiles:
    #TODO: get rid of this requirement
    if gff.file_type != 'gff':
      sys.exit('only gff supported: %s' % repr(gff))

    for gffeature in gff:
      # add the feature to the dict for remapping later
      flist = features.setdefault(gffeature.chrom, [])
      flist.append(gffeature)
      
      # consider only exons with Target attribute
      if gffeature.fields[2] != 'exon' or 'Target' not in gffeature.attrs: continue

      # extract target info  from the attributes
      target_id, target_start, target_end, target_strand = gffeature.attrs['Target'].split()
      
      # get or create locations dict for current read
      locs = readloc.setdefault(gffeature.chrom, dict())
      # get or create candidate location on target chromosome
      pos = locs.setdefault(target_id, target_start)
      if target_start < pos[target_id]:
        pos[target_id] = target_start

  # construct inverted map-list {chr -> [(location, read_name)]}
  invmap = dict()
  ambiguous = set()
  for read_name, locations in readloc.iteritems():
    if len(locations) > 1:
      ambiguous.add(read_name)
    else:
      seq, pos = locations.iteritems().next()
      read_list = invmap.setdefault(seq, [])
      read_list.append( (pos, read_name) )

  # compare chromosome ids
  # first by extracted numbers
  # then by length of the name
  # chr1 < chr1A < chrUn < chrUnV < chrZ
  import re
  num = re.compile('[0-9]+')
  def chromocompare(x, y):
    # extract any number from the chromosome name
    xmatch = num.search(x)
    ymatch = num.search(y)
    
    # case with one or no chromosome number in the comparison
    if not all([xmatch, ymatch]):
      if any([xmatch, ymatch]):
        return -cmp(xmatch, ymatch)
      else:
        return cmp(x, y)
    
    xnum = int(xmatch.group())
    ynum = int(ymatch.group())
    if cmp(xnum, ynum) == 0:
      return cmp(len(x), len(y))
    else:
      return cmp(xnum, ynum)
    

  chromosomes = invmap.keys()
  chromosomes.sort(cmp=chromocompare)

  # ambiguously assigned reads, assign position 0 for each of them
  invmap['chrAmb'] = [(0, read_name) for read_name in ambiguous]
  chromosomes.append('chrAmb']
  
  # unassigned reads
  all_reads = set(reads.iterkeys())
  
  

  # for the sorted list of chromosomes sort the positions and output the reads
  #FIXME: this stores all the sequences in memory (twice actually)
  # should not be a problem for cDNA 'genomes'
  vgenome = []
  for chr in chromosomes:
    # sort the positions
    invmap[chr].sort()
    chromseq = ""
    
    for pos, read_name in invmap[chr]:
      chromseq += reads[read_name]
      chromseq += 50 * 'N'
    
    vgenome.append((chr, chromseq))

if __name__ == "__main__": main()
