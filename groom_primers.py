#!/usr/bin/env python

# given the
# - annotations (gff, has to contain mRNA entries)
# - filtered variants (vcf, variants with PASS add to score of the locus)
# - primers in gff designed by design_primers.py
# produces 
# - cleaned up version of the gff with designed primers
# algorithm
# - go through the annotations
# - for each mRNA get all the primers/products
# - deduplicate all the entries based on position/strand (check primer3 values?)
# - apply the new naming scheme 
#   LU%04d[-%d(FRG)[-%d]] % contig_num, pcr_num, gt_num
#   where contig_num is nuber of contig in annotation
#         pcr_num is number of pcr product in the contig
#         gt_num is number of variant in pcr product
# - drop all genotyping primers with PROBLEMS from primer3
# - assign score to the pcr products based on number of PASS variants in the 
#   current mRNA (this should be a proxy for the level of site segregation)
# 
# Author: Libor Morkovsky, 2012
#

import sys
import os
import itertools
import string

import pysam
import vcf
import pybedtools

import operator
def deduplicated(l, _cmp=operator.eq):
    """return list with elements from l,
    but only with the last of the duplicit elements
    using given comparison function
    performs O(n^2) comparisons of elements
    """
    # test each of the combinations only once
    combinations = itertools.combinations(xrange(len(l)), 2)
    # go element by element
    c_by_elem = itertools.groupby(combinations, key=lambda e:e[0])
    # construct new list with elements not having any duplications
    result = [l[elem] for elem, it in c_by_elem if not any(_cmp(l[a], l[b]) for a, b in it)]
    # the last element is always unique
    result.append(l[-1])
    return result

def main():
    if len(sys.argv) < 4:
        sys.exit('use: %s mrna_gff_with_tbi variants_vcf_with_tbi primer_gff_with_tbi [starting_locus_id]' % sys.argv[0])

    annotations = pybedtools.BedTool(sys.argv[1])
    variants = vcf.Reader(filename=sys.argv[2])
    primers = pybedtools.BedTool(sys.argv[3])
    
    locus_ordinal = 0 if len(sys.argv == 4) else int(sys.argv[4])

    # go through all mRNA entries in the annotations
    # the other way around is to go variant by variant, and then filter out the already seen 
    # mrnas, which requires keeping a dict() of the seen mrnas
    for mrna in itertools.ifilter(lambda f: f.fields[2] == 'mRNA', annotations):
        # separate distinct types of features
        # (rejected-var, pcr-product, primer-pcr, primer-gt)
        mrna_primers = {}
        for f in primers.tabix_intervals(mrna):
            type_list = mrna_primers.setdefault(f.fields[2], [])
            type_list.append(f)

        # deduplicate each of the groups
        # if the position of the primers is equal, the sequences are 
        # also equal, so the termodynamic parameters should be the same as well
        # so we can use the Interval() '==' operator, which ignores gff attributes
        # the == operator raises error when the intervals contain each other
        def intrevals_eq(a, b):
            try:
                return a == b
            # when intervals contain each other, they're not equal
            except NotImplementedError:
                return False
                
        # skip the rejected vars
        mrna_primers_dd = {k:deduplicated(v, intrevals_eq) for k, v in mrna_primers.iteritems() if k != 'rejected-var'}

        # no primers in this mrna, fetch next one
        if len(mrna_primers_dd) == 0: 
            continue

        # count the interesting variants for current locus
        mrna_variants = variants.fetch(mrna.chrom, mrna.start, mrna.end)
        pass_variants = sum(1 for var in mrna_variants if not var.FILTER)

        # found contig with primers, increment locus number
        locus_ordinal += 1
        
        # decorate each pcr-product with the number of chosen variants 
        # in the locus
        for prod in mrna_primers_dd['pcr-product']:
            prod.attrs['variants_in_locus'] = pass_variants
        
        # remove the bad genotyping primers
        mrna_primers_dd['primer-gt'] = [p for p in mrna_primers_dd['primer-gt'] if 'PROBLEMS' not in p.attrs]
        

if __name__ == "__main__": main()
