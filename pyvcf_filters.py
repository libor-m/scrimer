#
# implementation of vcf filters
# for pyvcf vcf_filter.py
#
# Libor Morkovsky 2012
#
import vcf

class DepthPerSample(vcf.Filter):

  description = 'Threshold read depth per sample'
  name = 'dps'

  @classmethod
  def customize_parser(self, parser):
    parser.add_argument('--depth-per-sample', type=int, default=5,
              help='Filter out calls having lower per sample coverage [5]')

  def __init__(self, args):
    self.threshold = args.depth_per_sample

  def __call__(self, record):
    mindepth =  min([sam['DP'] for sam in record.samples])
    if mindepth < self.threshold:
      return mindepth
