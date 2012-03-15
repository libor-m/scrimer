#! /usr/bin/env python

#
# input is a gff file produced by the exon predictor software
#  and a set of bed files containing the features to be transfered
#  to the cDNA
# 
# output is an 'inverted' gff file describing the exons in the transcripts
#  and containing also all the features from the other input files, 
#  if they were overlapping with the predicted exons
#
# Author: Libor Morkovsky, 2012
#

import sys
import pybedtools
import itertools
import re

if len(sys.argv) < 2:
  sys.exit('use: %s gff_with_exons [bed_with_fetures] [...]' % sys.argv[0])
  
# load all feature files as IntervalFiles
# convert to bed6 if needed
bedtools = [pybedtools.BedTool(fname) for fname in sys.argv[2:]]
featurefiles = [bedtool.as_intervalfile() if bedtool.field_count() < 12 else bedtool.bed6().as_intervalfile() for bedtool in bedtools]

gff = pybedtools.BedTool(sys.argv[1])
if gff.file_type != 'gff':
  sys.exit('supply gff as the exon map file!')

#
# gmap output (GFF type), with annotations (GFF attribute):
#  gene - target id(Name)
#  mRNA - target id(Name), coverage (Coverage), identity (Identity)
#  exon - target position (Target)
#  CDS - target position (Target)
# sim4db output:
#  mRNA - target position (split(Name, ':')[1] and split(split(Target, ':')[1])[0]), targetLen, pA, pT, candidate geneRegion 
#  exon - target position (dtto), Gap, nMatches, intron (bad GFF fomat, can contain unescaped =)
# 
# --> to map gene and mRNA blocks from both formats (specifically gmap)
# we'd have either to fully parse the gff (with the relationships) - better solution
# or to get a faidx file for the target 'genome' (set of transcripts) 
# and set the mrna and gene features to contain the whole transcript - easier solution
# 
# go through the gff file, outputting the 'inverted' exons
# and potentially intersecting features from the other files
#
"""
GFF
  1: "seqid"
  2: "source"
  3: "type"
  4 & 5: "start" and "end"
  6: "score"
  7: "strand"
  8: "phase"
  9: "attributes"
"""
rmnum = re.compile('^[0-9]+:')

for gffeature in gff:
  # as a first approach we transfer only the exons and CDS
  if gffeature.fields[2] != 'exon' and gffeature.fields[2] != 'CDS': continue

  # extract target info  from the attributes
  target_id, target_start, target_end, target_strand = gffeature.attrs['Target'].split()
  # remove the leading #: if present (sim4db)
  target_id = rmnum.sub('', target_id)
  gffeature.chrom = rmnum.sub('', gffeature.chrom)
  
  # we do not transfer the parent features now, so do not break the format
  if 'Parent' in gffeature.attrs: del gffeature.attrs['Parent']
  
  # use the .fields[] to avoid coordinate conversion to bed format (-1 for start)
  gffeature.attrs['Target'] = "%s %s %s %s" % (gffeature.chrom, gffeature.fields[3], gffeature.fields[4], gffeature.strand)
  
  # transfer the source application from original gff to attributes in the output
  gffeature.attrs['source'] = gffeature.fields[1]
  
  # copy the rest of the attribues intact
  gflist = [target_id, 'liftover', gffeature.fields[2], 
    target_start, target_end, gffeature.score, target_strand,
    gffeature.fields[7], str(gffeature.attrs)]
  
  # create the inverted feature
  invgff = pybedtools.create_interval_from_list(gflist)
  print str(invgff).strip()

  # choose 
  if gffeature.fields[2] != 'exon': continue
  
  # find any features (consider them exons) in the additional files
  # that intersect with current exon (merge the hits from all files)
  intersecting = itertools.chain.from_iterable([ff.all_hits(gffeature) for ff in featurefiles])
  
  # for each intersecting feature create 
  for feature in intersecting: