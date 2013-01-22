"""
Implementation of vcf filters for pyvcf `vcf_filter.py`.

Author: Libor Morkovsky 2012
"""

# This file is a part of Scrimer.
# See LICENSE.txt for details on licensing.
#    Copyright (C) 2012, 2013 Libor Morkovsky

import vcf

class DepthPerSample(vcf.Filter):

  description = 'Threshold read depth per sample'
  name = 'dps'

  @classmethod
  def customize_parser(self, parser):
    parser.add_argument('--depth-per-sample', type=int, default=5,
              help='Minimum required coverage in each sample')

  def __init__(self, args):
    self.threshold = args.depth_per_sample

  def __call__(self, record):
    # do not test depth for indels
    if record.is_indel: return
    
    mindepth =  min([sam['DP'] for sam in record.samples])
    if mindepth < self.threshold:
      return mindepth

class AvgDepthPerSample(vcf.Filter):

  description = 'Threshold average read depth per sample (read_depth / sample_count)'
  name = 'avg-dps'

  @classmethod
  def customize_parser(self, parser):
    parser.add_argument('--avg-depth-per-sample', type=int, default=3,
              help='Minimum required average coverage per sample')

  def __init__(self, args):
    self.threshold = args.avg_depth_per_sample

  def __call__(self, record):
    avgcov = float(record.INFO['DP']) / len(record.samples)
    if avgcov < self.threshold:
      return avgcov

class SnpOnly(vcf.Filter):
  
  description = 'Choose only SNP variants'
  name = 'snp-only'

  def __call__(self, record):
    if not record.is_snp: return True

  def filter_name(self):
      return self.name

class DistinguishingVariants(vcf.Filter):

  description = 'Given a group of samples, choose variants that are not shared with the rest of the samples'
  name = 'contrast-samples'

  @classmethod
  def customize_parser(self, parser):
    #TODO: argparse.add_mutually_exclusive_group
    # for --sample-numbers
    parser.add_argument('--sample-names', nargs='+', metavar='sample_id',
              help='Names of samples in the group to contrast against the rest')

  def filter_name(self):
      return self.name

  def __init__(self, args):
    self.group = set(args.sample_names)
    import re
    self.splitter = re.compile('[/|]')


  def __call__(self, record):
    # find the names for rest of the samples
    #TODO: this could be done only once in __init__, if there was access to the reader
    rest = set(sam.sample for sam in record.samples if sam.sample not in self.group)
    
    # create sets of alleles in group and in the rest
    from itertools import chain
    alleles = lambda grp: set(chain.from_iterable(self.splitter.split(record.genotype(sam).gt_bases) for sam in grp))
    
    # if there are any alleles in common, filter the record
    common = alleles(rest) & alleles(self.group)
    if len(common): return common

