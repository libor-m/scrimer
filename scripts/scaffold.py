#! /usr/bin/env python

"""
Input

- fasta file with the original sequences
- set of gff files with exon features having the Target attribute 
  (product of ``liftover.py``)
- other gff/bed files to remap to the genome

Output 

- fasta with the input sequences pasted with 100 Ns

  - unambiguos sequences assigned to 'chromosomes' in the order 
    of the template genome (that was used to generate the Target exon mappings)
  - ambiguos (having more than one possible chromosome) assigned to **chrAmb**
  - unmapped sequences assigned to **chrUnmapped** (to distinguish from chrUn in target genome)

- gff file with the locations of the input sequences (gene)
  and remapped contents of the input gff files

Algorithm

- go through the input gff files, construct a dictionary {read_name -> {chr -> location}}
- add the lowest coordinate found on a given chromosome (overwriting previous values)
- sort the list with single candidate locations with 'chromocompare'

Author: Libor Morkovsky, 2012
"""

# This file is a part of Scrimer.
# See LICENSE.txt for details on licensing.
#    Copyright (C) 2012, 2013 Libor Morkovsky

import sys
import pybedtools

# reuse the generator fasta reader from lastz
from fasta_fragments import fasta_sequences

def get_overlap(a, b):
  return max(0, min(a[1], b[1]) - max(a[0], b[0]))
  
def merge_overlapping(ranges):
  """
  for a list of intervals stored as tuples
  return a list of intervals where 
  the overlapping regions are merged
  """
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
      # if there is an overlap or touch
      # try to extend current region
      if nstart <= end:
        end = max(end, nend)
        # continue with merging after this element
        idx = nidx
      # we're in sorted list, there will be no more overlaps
      else:
        break
    
    idx += 1
    newranges.append((start, end))
  
  return newranges
  
def coverage_by_target(contig, features, type="exon"):
  """
  for all features of given type
  sum the coverage with respect to the
  Target attribute
  return a dictionary {target_id: fraction_of_contig_covered}
  """
  # harvest the regions
  regions = dict()
  result = dict()
  for gffeature in features:
    if gffeature.fields[2] != 'exon' or 'Target' not in gffeature.attrs: continue

    # extract target info from the attributes
    target_id = gffeature.attrs['Target'].split()[0]
    
    # add to regions
    range_list = regions.setdefault(target_id, [])
    range_list.append( (int(gffeature.start), int(gffeature.end)) )
    
  # merge overlapping regions 
  for seq, ranges in regions.iteritems():
    newranges = merge_overlapping(ranges)
    
    # calculate the coverage
    covered_bases = sum(map(lambda x: get_overlap((contig.start, contig.end), x), newranges))
  
    result[seq] = float(covered_bases) / (contig.end - contig.start + 1)
    
  return result
    
  
def fasta_write(file, seqlist, line_len=70):
  for name, seq in seqlist:
    file.write(">%s\n" % name)
    for i in xrange(0, len(seq), line_len):
      file.write("%s\n" % seq[i : i + line_len])

import re
num = re.compile('[0-9]+')
def chromocompare(x, y):
  # compare chromosome ids
  # first by extracted numbers
  # then by length of the name
  # chr1 < chr1A < chrUn < chrUnV < chrZ

  # extract any number from the chromosome name
  xmatch = num.search(x)
  ymatch = num.search(y)
  
  # case with one or no chromosome number in the comparison
  if not all([xmatch, ymatch]):
    if any([xmatch, ymatch]):
      # None < anything, we need it >, return a negation of the comparison
      return -cmp(xmatch, ymatch)
    else:
      return cmp(x, y)
  
  xnum = int(xmatch.group())
  ynum = int(ymatch.group())
  if cmp(xnum, ynum) == 0:
    return cmp(len(x), len(y))
  else:
    return cmp(xnum, ynum)

def main():
  if len(sys.argv) < 2:
    sys.exit('use: %s contigs.fasta gff_with_mapped_exons [gff_with_mapped_exons ...] output.fasta remapped.gff' % sys.argv[0])
  
  # load the sequences into dict {read_name -> sequence}
  print >> sys.stderr, "loading sequences from %s" % sys.argv[1]
  seqfile = open(sys.argv[1])
  reads = { read_name : seq for (read_name, seq) in fasta_sequences(seqfile) }
  seqfile.close()
  
  # construct {read_name -> {chr -> location}}
  readloc = dict()
  features = dict()
  gffiles = [pybedtools.BedTool(fname) for fname in sys.argv[2:-2]]
  for gff in gffiles:
    #TODO: get rid of this requirement
    if gff.file_type != 'gff':
      sys.exit('only gff supported: %s' % repr(gff))

    print >> sys.stderr, "parsing gff %s" % repr(gff)

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
      pos = locs.setdefault(target_id, int(target_start))
      if target_start < pos:
        locs[target_id] = int(target_start)

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
    
  chromosomes = invmap.keys()
  chromosomes.sort(cmp=chromocompare)

  # try to resolve the ambiguity by exon coverage
  print >> sys.stderr, "trying to resolve %d ambiguities " % len(ambiguous),
  amb_resolved = set()
  for read_name in ambiguous:
    coverages = [(coverage, target) for target, coverage in coverage_by_target(pybedtools.Interval(read_name, 0, len(reads[read_name])), features[read_name]).iteritems()]
    coverages.sort(reverse=True)

    # get an idea about the disambiguation process
    # print "disamb:", read_name, coverages
    
    # do the decisions only for decently covered contigs, leave others as amb
    if coverages[0][0] < 0.5: continue
    
    # if the best hit covers 4 times more than the second one, use it
    if coverages[0][0] < 4 * coverages[1][0] : continue

    # append the read to the chosen chromosome, 
    # as the others
    pos = readloc[read_name][coverages[0][1]]
    read_list = invmap.setdefault(coverages[0][1], [])
    read_list.append( (pos, read_name) )
    
    amb_resolved.add(read_name)

  print >> sys.stderr, "%d resolved" % len(amb_resolved)
  
  # ambiguously assigned reads, assign position 0 to all
  ambiguous -= amb_resolved
  invmap['chrAmb'] = [(0, read_name) for read_name in ambiguous]
  chromosomes.append('chrAmb')
  
  # unassigned reads
  all_reads = set(reads.iterkeys())
  assigned_reads = set(read_name for chrom, read_list in invmap.iteritems() for pos, read_name in read_list)
  unassigned_reads = all_reads - assigned_reads
  invmap['chrUnmapped'] = [(0, read_name) for read_name in unassigned_reads]
  chromosomes.append('chrUnmapped')

  # for the sorted list of chromosomes sort the positions and output the reads
  print >> sys.stderr, "writing output to %s and %s" % (sys.argv[-2], sys.argv[-1])
  splitlen = 150
  splitter = splitlen * 'N'

  fo_gff = open(sys.argv[-1], "w")
  fo_fasta = open(sys.argv[-2], "w")

  for chr in chromosomes:
    # sort the positions
    invmap[chr].sort()
    
    # construct the chromosome sequence
    chromseq = splitter.join(reads[read_name] for pos, read_name in invmap[chr])

    # output it as fasta
    fasta_write(fo_fasta, [(chr, chromseq)])
    
    # output a gff mRNA feature for each contig put on the chromosome
    # and remap all gff features from the input for reads at current chr
    feature_lengths = map(len, (reads[read_name] for pos, read_name in invmap[chr]))
    feature_positions = map(lambda x: sum(feature_lengths[:x]) + splitlen * x, xrange(0, len(feature_lengths)))
    for (index, (pos, read_name)) in enumerate(invmap[chr]):
      
      at = pybedtools.Attributes(' ');
      at['ID'] = read_name
      flist = [chr, 'virtual_genome', 'mRNA',
               str(feature_positions[index] + 1), str(feature_positions[index] + feature_lengths[index]),
               '1', '+', '.', str(at)]
      
      fo_gff.write(str(pybedtools.create_interval_from_list(flist)))
      
      # remap all the other features
      if read_name not in features: continue
      for feature in features[read_name]:
        # feature.attrs['Parent'] = read_name
        feature.chrom = chr
        feature.start += feature_positions[index]
        feature.end += feature_positions[index]
        fo_gff.write(str(feature))

  fo_gff.close()
  fo_fasta.close()
 
  
if __name__ == "__main__": main()
