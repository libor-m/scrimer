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
              help='Minimum required coverage in each sample [5]')

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
              help='Minimum required average coverage per sample [3]')

  def __init__(self, args):
    self.threshold = args.avg_depth_per_sample

  def __call__(self, record):
    avgcov = float(record.INFO['DP']) / len(record.samples)
    if avgcov < self.threshold:
      return avgcov
