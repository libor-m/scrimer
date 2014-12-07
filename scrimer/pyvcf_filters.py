"""
Implementation of vcf filters for :ref:`pyvcf <PyVCF>` `vcf_filter.py`.

Author: Libor Morkovsky 2012
"""

# This file is a part of Scrimer.
# See LICENSE.txt for details on licensing.
#    Copyright (C) 2012, 2013 Libor Morkovsky

import vcf

class DistinguishingVariants(vcf.Filter):
    """Given a group of samples, choose variants that 
    are not shared with the rest of the samples
    """
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
      try:
          common = alleles(rest) & alleles(self.group)
      except:
          return 'contrast-exception'

      if len(common): return common

