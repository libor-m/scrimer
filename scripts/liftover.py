#! /usr/bin/env python

"""
Input 

- gff file produced by the exon predictor software
- set of bed files containing the features to be transfered
  to the cDNA

Output

- 'inverted' gff file annotating exons in the transcripts
  and containing also all the features from the other input files, 
  if they were overlapping with the predicted exons

Author: Libor Morkovsky, 2012
"""

# This file is a part of Scrimer.
# See LICENSE.txt for details on licensing.
#    Copyright (C) 2012, 2013 Libor Morkovsky

import sys
import itertools
import re

import pybedtools

def main():
    if len(sys.argv) < 2:
      sys.exit('use: %s gff_with_exons [bed_with_fetures] [...]' % sys.argv[0])
      
    # load all feature files as IntervalFiles
    # convert to bed6 if needed
    bedtools = [pybedtools.BedTool(fname) for fname in sys.argv[2:]]
    # fix for some damned bug in pybedtools
    for bt in bedtools:
        bt._isbam = False
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
      
      # sanity check, do not go on if coords are strange
      if int(target_start) > int(target_end): 
        print >> sys.stderr, "target coords weird (%d, %d)" % (int(target_start), int(target_end))
        continue
      
      # remove the leading 'number:' if present (sim4db)
      target_id = rmnum.sub('', target_id)
      gffeature.chrom = rmnum.sub('', gffeature.chrom)
      
      # we do not transfer the parent features now, 
      # so do not break the format and remove the Parent entry
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

      # look for annotations overlapping only exons (CDS are contained in exons) 
      if gffeature.fields[2] != 'exon': continue
      
      # find features (consider them exons) in the additional files
      # that intersect with current exon (merge the hits from all files)
      intersecting = itertools.chain.from_iterable([ff.all_hits(gffeature) for ff in featurefiles])
      
      # for each intersecting feature create a clipped feature
      for feature in intersecting:
        
        fs, fe = feature.start, feature.end

        #FIXME: remove the ' ' arg, temporary fix for pybedtools
        at = pybedtools.Attributes(' ');
        
        
        at['Name'] = feature.name
        at['coords'] = gffeature.fields[1]
        
        # this is not necessary if the extension is .gff3 or there is a #version tag in the 
        # beginning of the gff file
        # -
        # because of the igv feature merging, add some uniqness 
        #at['seqid'] = target_id
        
        # left clip
        lclip = gffeature.start - feature.start
        if lclip > 0:
          at['leftClip'] = str(lclip)
          feature.start = int(target_start)
        else:
          feature.start = max(int(target_start) - lclip, 1)
        
        # right clip
        rclip = feature.end - gffeature.end
        if rclip > 0:
          at['rightClip'] = str(rclip)
          feature.end = int(target_end)
        else:
          feature.end = min(int(target_start) + (feature.end - gffeature.start), int(target_end))
        
        #TODO: review the rebasing calculation, still produces start>end
        # quick fix:
        feature.start, feature.end = sorted([feature.start, feature.end])
        
        # sanity check
        # if feature.start < 1:
        #  print >> sys.stderr, str(feature).strip(), target_start, target_end, lclip, rclip, fs, fe
        
        flist = [target_id, 'liftover', 'exon',
          str(feature.start), str(feature.end), str(gffeature.score), feature.strand,
          '.', str(at)]
        
        print str(pybedtools.create_interval_from_list(flist)).strip()

if __name__ == "__main__": main()
