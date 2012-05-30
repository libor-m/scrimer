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
# - deduplicate all the entries based on position/strand (and check primer3 values)
# - apply the new naming scheme 
#   LU%04d[-%d(FRG)[-%d]] % contig_num, pcr_num, gt_num
#   where contig_num is nuber of contig in annotation
#         pcr_num is number of pcr product in the contig
#         gt_num is number of variant in pcr product
# - drop all color=#bb0000 entries (errors)
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
    if len(sys.argv) < 2:
        sys.exit('use: %s mrna_gff_with_tbi variants_vcf_with_tbi primer_gff_with_tbi' % sys.argv[0])

    annotations = pybedtools.BedTool(sys.argv[1])
    variants = vcf.Reader(filename=sys.argv[2])
    primers = pybedtools.BedTool(sys.argv[3])
    # open the variants once more, so the ifilter is not broken by .fetch

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

            
        # get a minimal predicted exon
        # (trying to be conservative)
        # consider only predicted, not transferred exons (source in attrs)
        ex_predicted = [f for f in var_exons if 'source' in f.attrs]
        ex_start = max(f.start for f in ex_predicted)
        ex_end = min(f.end for f in ex_predicted)
        min_exon = pybedtools.Interval(var.CHROM, ex_start, ex_end)
        
        # check if the exon has a sane size (> 70)
        if min_exon.length < 70: 
            report_rejected_var('minimal exon too short for PCR')
            continue

        # check if the variant is in position that permits PCR amplification
        # (more than 20 bases from both ends)
        # this implies a check if the variant is still inside of the minimal exon
        distances = (var.POS - min_exon.start, min_exon.end - var.POS)
        if min(distances) < 20: 
            report_rejected_var('variant too close to boundary or outside of minimal exon')
            continue
        
        # get all variants for the minimal exon
        exon_variants = var_fetcher.fetch(var.CHROM, min_exon.start, min_exon.end)
        
        # get sequence for the minimal exon
        # and patch all the variants with Ns
        seq = genome.fetch(min_exon.chrom, min_exon.start, min_exon.end)
        pseq = patched_sequence(seq, min_exon.start, exon_variants)
        
        # check if there is at least 20 fixed bases on either side
        # of the target variant because of the genotyping primer
        # (this could be faster exit path than a call to primer3)
        gt_primer_seqs = get_flanking(pseq, var.POS - 1 - min_exon.start, max_gt_primer_len)
        # attach 'N' to end of each primer, so we get the full len
        # instead of -1, when there is no N in the original primer
        max_gt_lens = [(s + 'N').find('N') for s in  [gt_primer_seqs[0][::-1], gt_primer_seqs[1]]]
        if all(x < min_gt_primer_len for x in max_gt_lens): 
            report_rejected_var('no possible genotyping primer longer than %d bp' % min_gt_primer_len)
            continue
        
        # call primer3 to find suitable genotyping primers
        gt_primers = find_gt_primers(pseq, var.POS - 1 - min_exon.start)
        if all(len(p) == 0 for p in [gt_primers[0]['LEFT'], gt_primers[0]['RIGHT']]):
            report_rejected_var('no good genotyping primer found by primer3')
            continue

        # call primer3 to design PCR primers
        # mark the region with suitable genotyping primers as a target
        pcr_min = pcr_max = var.POS - 1 - min_exon.start
        if len(gt_primers[0]['LEFT']):
            pcr_min = int(gt_primers[0]['LEFT'][0]['position'].split(',')[0])
        if len(gt_primers[0]['RIGHT']):
            # primer3 uses coordinates of the 5' end
            pcr_max = int(gt_primers[0]['RIGHT'][0]['position'].split(',')[0])
            
        #TODO: call primer3 only once with multiple records?
        primers = find_pcr_primers(pseq, (pcr_min, pcr_max - pcr_min))
        
        # generate the name for the primer ensemble
        # use the name from reference genome only if it is unique
        # otherwise use the transcript name
        ref_names = set(e.attrs['Name'] for e in var_exons if 'coords' in e.attrs)
        if len(ref_names) == 1:
            primer_name = ref_names.pop()
        else:
            names = [e.attrs['Name'] for e in ex_predicted if 'Name' in e.attrs]
            primer_name = names[0] if len(names) else 'NAME-UNKNOWN'
        
        name_ordinal = primer_names.setdefault(primer_name, 1)
        primer_names[primer_name] += 1
        
        # if something was found, output the first pair
        # (should be the best one)
        if len(primers[0]['PAIR']):
            prod_id = 'LU-%s-%d' % (primer_name, name_ordinal)
            # calculate the 'position' entry PAIR, so it is the same as in LEFT/RIGHT entries and can be used by primer_to_gff
            # to represent the whole-product
            primers[0]['PAIR'][0]['position'] = primers[0]['LEFT'][0]['position'].split(',')[0] + ',' + primers[0]['PAIR'][0]['PRODUCT_SIZE']
            print primer_to_gff(prod_id, primers[0]['PAIR'][0], 'pcr-product', var.CHROM, min_exon.start, '+', locus_name=primer_name)
            print primer_to_gff('LU-PCR-%s-%dF' % (primer_name, name_ordinal), primers[0]['LEFT'][0],  'primer-pcr', var.CHROM, min_exon.start, '+', Parent=prod_id)
            print primer_to_gff('LU-PCR-%s-%dR' % (primer_name, name_ordinal), primers[0]['RIGHT'][0], 'primer-pcr', var.CHROM, min_exon.start, '-', Parent=prod_id)
        else:
            report_rejected_var('no suitable PCR primers found')
            continue
        
        if len(gt_primers[0]['LEFT']):
            color = '#bb0000' if 'PROBLEMS' in gt_primers[0]['LEFT'][0] else '#00bb00'
            print primer_to_gff('LU-GT-%s-%dF' % (primer_name, name_ordinal), gt_primers[0]['LEFT'][0], 'primer-gt', var.CHROM, min_exon.start, '+', color=color, Parent=prod_id)
        if len(gt_primers[0]['RIGHT']):
            color = '#bb0000' if 'PROBLEMS' in gt_primers[0]['RIGHT'][0] else '#00bb00'
            print primer_to_gff('LU-GT-%s-%dR' % (primer_name, name_ordinal), gt_primers[0]['RIGHT'][0], 'primer-gt', var.CHROM, min_exon.start, '-', color=color, Parent=prod_id)
        
        #TODO: here should be some kind of scoring for the primer ensemble
        # - number of other segregating polymorphisms in the same exon/?mRNA 
        # - variant overall coverage, minimal sample coverage - from vcf
        # - 

if __name__ == "__main__": main()
