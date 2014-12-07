#!/usr/bin/env python

"""
Input

- reference sequence (fasta with samtools fai index)
- annotations (gff3, has to contain exon entries)
- filtered variants (vcf, primers are designed for variants with PASS)

Output

- PCR and genotyping primers selected using primer3 (gff3)

Algorithm

- there is only a few selected variants, so the least amount of work
  will be to do the work only for variants
- for each of the selected variants

  - request exons
  - apply the technical constraints 
    (minimal primer length of 20 from the edge of an exon)
  - patch exon sequence to mark positions of known variants
  - find suitable genotyping primers
  - design PCR primers to flank the (usable) genotyping primers

Author: Libor Morkovsky, 2012, 2014
"""

# This file is a part of Scrimer.
# See LICENSE.txt for details on licensing.
#    Copyright (C) 2012, 2013 Libor Morkovsky

import sys
import os
import itertools
import argparse

import pysam
import vcf
import pybedtools

from scrimer import primer3_connector

min_gt_primer_len = 20
pref_gt_primer_len = 22
max_gt_primer_len = 28

# get location of config dir from environment
if 'PRIMER3_CONFIG' in os.environ:
    PRIMER3_CONFIG = os.environ['PRIMER3_CONFIG']
else:
    PRIMER3_CONFIG = '/opt/primer3/bin/primer3_config/'

def patched_sequence(seq, seq_start, variants):
    """Returns sequence with all variant sites replaced by Ns
    a functional approach forced by python immutable strings
    seq - sequence string
    seq_start - start of sequence in refernce for variants
    variants - [vcf._Record]
    """
    # consider each variable site only once
    uniq_positions = set((var.POS - 1 - seq_start) for var in variants)
    
    # break the input seq apart, leaving out the variant sites
    ends = sorted(uniq_positions)
    begins = [0] + [pos + 1 for pos in ends]
    ends.append(len(seq))
    fragments = (seq[b:e] for (b, e) in zip(begins, ends))
    
    # substite the variant sites with Ns
    return 'N'.join(fragments)

#TODO: consider PRIMER_LOWERCASE_MASKING=1
def find_pcr_primers(pseq, target):
    """call primer3 executable, return the results
    """
    def_params = {
    'PRIMER_THERMODYNAMIC_PARAMETERS_PATH': PRIMER3_CONFIG,
    'PRIMER_MAX_NS_ACCEPTED':'0',
    }
    
    p3 = primer3_connector.Primer3(**def_params)

    rec_params = {
    'SEQUENCE_ID':'pick_pcr_primers_task',
    'SEQUENCE_TEMPLATE':pseq,
    'SEQUENCE_TARGET': "%d,%d" % (target[0], target[1]),
    'PRIMER_OPT_SIZE':'19',
    'PRIMER_MIN_SIZE':'17',
    'PRIMER_MAX_SIZE':'25',
    'PRIMER_PRODUCT_SIZE_RANGE':'70-300',
    'PRIMER_LOWERCASE_MASKING':'1',
    }
    
    return p3.call([rec_params])

def find_gt_primers(pseq, target):
    """call primer3 executable, return the results
    simulating pick_discriminating_primers with FORCE_*_END
    to avoid off-by-one error present in primer3
    """
    def_params = {
    'PRIMER_THERMODYNAMIC_PARAMETERS_PATH': PRIMER3_CONFIG,
    'PRIMER_MAX_NS_ACCEPTED':'0',
    }
    
    p3 = primer3_connector.Primer3(**def_params)

    rec_params = {
    'PRIMER_TASK':'generic',
    'SEQUENCE_ID':'pick_gt_primers_task',
    'SEQUENCE_TEMPLATE': pseq,
    'SEQUENCE_TARGET': "%d,1" % target,
    'SEQUENCE_FORCE_LEFT_END': str(target - 1),
    'SEQUENCE_FORCE_RIGHT_END': str(target + 1),
    'PRIMER_OPT_SIZE': str(pref_gt_primer_len),
    'PRIMER_MIN_SIZE': str(min_gt_primer_len),
    'PRIMER_MAX_SIZE': str(max_gt_primer_len),
    'PRIMER_PRODUCT_SIZE_RANGE':'%d-%d' % (2*min_gt_primer_len + 1, 2*max_gt_primer_len + 1),
    'PRIMER_MIN_TM':'50',
    'PRIMER_OPT_TM':'55',
    'PRIMER_NUM_RETURN':'1',
    # pick the primers even if there is no decent pair, we don't need pairs
    'PRIMER_PICK_ANYWAY': '1', 
    # pick a thermodynamically worse primer, but avoid Ns at all costs
    'PRIMER_WT_NUM_NS': '100', 
    }
    
    return p3.call([rec_params])

def primer_to_gff(name, primer, tag, seq_name, seq_start, strand, **kwargs):
    """Create a gff feature from 
    partially parsed primer3 results.
    """
    pos, len = map(int, primer['position'].split(','))
    
    # transfer the calculated values to attributes
    # skip the fields used elsewhere in gff
    at = pybedtools.Attributes(' ')
    for k, v in primer.iteritems():
        if k == 'position': continue
        at[k] = v.replace(';', '%3B')
        
    at['ID'] = name
    
    # pass all optional params to attributes
    at.update({k:str(v) for k, v in kwargs.iteritems()})
    
    # primer3 provides the coordinates of right primer with the 
    # pos pointing to the last base
    if strand == '-':
        start = seq_start + pos - len + 2
        end = seq_start + pos + 1
    else:
        start = seq_start + pos + 1
        end = seq_start + pos + len

    gflist = [seq_name, 'design-primers', tag, 
        str(start), str(end),
        primer['PENALTY'], strand, '.', str(at)]

    return pybedtools.create_interval_from_list(gflist)

def get_flanking(seq, pos, size):
    """Get sequence flanking given position
    add Ns if there is no input data (ends of input sequence)
    """
    left = seq[max(pos - size, 0):pos]
    right = seq[pos + 1:pos + size + 1]
    return (left.rjust(size, 'N'), right.ljust(size, 'N'))

def reverse_complement(seq):
    """Return a DNA reverse complement
    """
    compl = string.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(compl)[::-1]

def primers_for_var(genome, annotations, variants_random_access, var):
    """In ideal case this returns list of gff features representing
    - a pair of pcr primers
    - pcr product (carries the thermodynamic penalties for using 
        the combination of current pcr primers)
    - two genotyping primers
    or a feature carrying a Note about the error
    """
    
    def rejected_var_feature(reason):
        """ Output a dummy feature holding the 
        reason for rejection of selected variant in 
        the design process.
        """
        at = pybedtools.Attributes(' ')
        at['Note'] = reason
        at['color'] = "#bb0000"
        
        gflist = [var.CHROM, 'design-primers', 'rejected-var', 
            str(var.POS - max_gt_primer_len), str(var.POS + max_gt_primer_len),
            '0', '+', '.', str(at)]
            
        return pybedtools.create_interval_from_list(gflist)

    # container for the output
    result = []

    var_features = annotations.tabix_intervals(pybedtools.Interval(var.CHROM, var.POS, var.POS + 1))
    var_exons = [f for f in var_features if f.fields[2] == 'exon']
    

    if len(var_exons) == 0: 
        return [rejected_var_feature('no predicted exons')]
    
    # get a minimal predicted exon
    # (trying to be conservative)
    # consider only predicted, not transferred exons (source in attrs)
    ex_predicted = [f for f in var_exons if 'source' in f.attrs]
    ex_start = max(f.start for f in ex_predicted)
    ex_end = min(f.end for f in ex_predicted)
    min_exon = pybedtools.Interval(var.CHROM, ex_start, ex_end)
    
    # check if the exon has a sane size (> 70)
    if min_exon.length < 70: 
        return [rejected_var_feature('minimal exon too short for PCR')]

    # check if the variant is in position that permits PCR amplification
    # (more than 20 bases from both ends)
    # this implies a check if the variant is still inside of the minimal exon
    distances = (var.POS - min_exon.start, min_exon.end - var.POS)
    if min(distances) < 20: 
        return [rejected_var_feature('variant too close to boundary or outside of minimal exon')]
    
    # get all variants for the minimal exon
    exon_variants = variants_random_access.fetch(var.CHROM, min_exon.start, min_exon.end)
    
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
        return [rejected_var_feature('no possible genotyping primer longer than %d bp' % min_gt_primer_len)]

    
    # call primer3 to find suitable genotyping primers
    gt_primers = find_gt_primers(pseq, var.POS - 1 - min_exon.start)
    if all(len(p) == 0 for p in [gt_primers[0]['LEFT'], gt_primers[0]['RIGHT']]):
        return [rejected_var_feature('no good genotyping primer found by primer3')]

    # call primer3 to design PCR primers
    # mark the region with suitable genotyping primers as a target
    pcr_min = pcr_max = var.POS - 1 - min_exon.start
    if len(gt_primers[0]['LEFT']):
        pcr_min = int(gt_primers[0]['LEFT'][0]['position'].split(',')[0])
    if len(gt_primers[0]['RIGHT']):
        # primer3 uses coordinates of the 5' end
        pcr_max = int(gt_primers[0]['RIGHT'][0]['position'].split(',')[0])
        
    
    # find primers in the sequence patched with Ns
    primers = find_pcr_primers(pseq, (pcr_min, pcr_max - pcr_min))
    
    # if something was found, output the first pair
    # (should be the best one)
    if len(primers[0]['PAIR']):
        # calculate the 'position' entry PAIR, so it is the same as in LEFT/RIGHT entries and can be used by primer_to_gff
        # to represent the whole-product
        primers[0]['PAIR'][0]['position'] = primers[0]['LEFT'][0]['position'].split(',')[0] + ',' + primers[0]['PAIR'][0]['PRODUCT_SIZE']
        result.append(primer_to_gff('pcr-product', primers[0]['PAIR'][0], 'pcr-product', var.CHROM, min_exon.start, '+'))
        result.append(primer_to_gff('pcr-left', primers[0]['LEFT'][0],  'pcr-primer', var.CHROM, min_exon.start, '+'))
        result.append(primer_to_gff('pcr-right', primers[0]['RIGHT'][0], 'pcr-primer', var.CHROM, min_exon.start, '-'))
    else:
        return [rejected_var_feature('no suitable PCR primers found')]
    
    # decorate all genotyping primers with basic statistics about the variant being genotyped
    
    vkeys = ['FQ', 'MQ']
    more = {('VAR_%s' % k):var.INFO[k] for k in vkeys if k in var.INFO}
    more['VAR_mindps'] = min(sam['DP'] for sam in var.samples)
    if 'DP4' in var.INFO:
        more['VAR_dp4'] = ','.join(map(str, var.INFO['DP4']))

    if len(gt_primers[0]['LEFT']):
        color = '#bb0000' if 'PROBLEMS' in gt_primers[0]['LEFT'][0] else '#00bb00'
        result.append(primer_to_gff('gt-left', gt_primers[0]['LEFT'][0], 'gt-primer', var.CHROM, min_exon.start, '+', 
            color=color, **more))
    if len(gt_primers[0]['RIGHT']):
        color = '#bb0000' if 'PROBLEMS' in gt_primers[0]['RIGHT'][0] else '#00bb00'
        more = {('VAR_%s' % k):var.INFO[k] for k in vkeys if k in var.INFO}
        result.append(primer_to_gff('gt-right', gt_primers[0]['RIGHT'][0], 'gt-primer', var.CHROM, min_exon.start, '-', 
            color=color, **more))
            
    return result

def argparser():
    args = argparse.ArgumentParser(description="Design primers using PRIMER3.")

    args.add_argument("--primer-min", type=int, default=20,
        help="Minimal length of the genotyping primer (default: 20).")
    args.add_argument("--primer-pref", type=int, default=22,
        help="Preferred length of the genotyping primer (default: 22).")
    args.add_argument("--primer-max", type=int, default=28,
        help="Minimal length of the genotyping primer (default: 28).")

    #args.add_argument("--annot-type", default="exon",
    #    help="Annotatoin type that should represent contiguous sequence in DNA (default: exon).")

    args.add_argument("genome",
        help="Fasta file with associated .fai index (use samtools fai to index).")
    args.add_argument("annots",
        help="GFF file with associated .tbi (use tabix to index).")
    args.add_argument("variants",
        help="VCF file with selected variants.")

    return args

def crc(line):
    return zlib.crc32(line)
    
def main():

    parser = argparser()
    args = parser.parse_args() 

    genome = pysam.Fastafile(args.genome)
    annotations = pybedtools.BedTool(args.annots)
    variants = vcf.Reader(filename=args.variants)
    # open the variants once more, so the ifilter is not broken by .fetch
    var_fetcher = vcf.Reader(filename=args.variants)

    # the most commonly changed value will be 
    # preferred length, fix the min/max if pref
    # is not in the default range
    if args.primer_pref < args.primer_min:
        args.primer_min = args.primer_pref
    if args.primer_pref > args.primer_max:
        args.primer_max = args.primer_pref
        
    min_gt_primer_len = args.primer_min
    pref_gt_primer_len = args.primer_pref
    max_gt_primer_len = args.primer_max

    # locus can be either exon or mrna, for now we consider mrna a locus
    # locus_feature_type = args.annot_type
    
    seen_loci = set()

    # by-locus naming convention implies following algo
    # (almost single pass, keeping only the names of loci already seen 
    #  in memory)
    # algorithm:
    # get a chosen variant
    # find locus feature for chosen variant
    # if locus already seen, get next variant
    # find all other chosen variants for locus
    # for each variant in locus design primers
    #  if designed primers overlap with others in locus, deduplicate
    
    # use iterator filter, so the variants are streamed
    for var in itertools.ifilter(lambda v: not v.FILTER, variants):
        # get a locus feature for the var
        var_features = annotations.tabix_intervals(pybedtools.Interval(var.CHROM, var.POS, var.POS + 1))
        # var_loci = [f for f in var_features if f.fields[2] == locus_feature_type]
        var_loci = [f for f in var_features if f.fields[2] == 'mRNA']
        
        # var is not in any known feature
        if len(var_loci) == 0:
            continue

        # quick check for now
        if len(var_loci) > 1:
            print >> sys.stderr, 'too many locus features for %s: %s' % (str(var), str(var_loci))
            continue
        
        locus = var_loci[0]
        
        # all variants from this locus are already done
        if locus.name in seen_loci:
            continue
        else:
            seen_loci.add(locus.name)

        # process all variants in current locus
        # set the names and parenting structure
        locus_variants = [v for v in var_fetcher.fetch(var.CHROM, locus.start, locus.end)]
        locus_primers = []
        pcr_id = 0
        locus_out_name = "LU%04d" % len(seen_loci)
        
        # summary statistics for whole locus
        # used to decorate all products
        pass_variants = sum(1 for var in locus_variants if not var.FILTER)
        
        for lvar_num, lvar in enumerate(itertools.ifilter(lambda v: not v.FILTER, locus_variants)):
            # create a convinient dict from the max-5 feature result
            curr_primers = {f.name:f for f in primers_for_var(genome, annotations, var_fetcher, lvar)}

            # if there was an error designing the primers
            # output the error indicating features
            if 'pcr-product' not in curr_primers:
                map(lambda f: sys.stdout.write(str(f)), curr_primers.itervalues())
                continue

            def add_feature(f, name, parent):
                f.name = name
                f.attrs['Parent'] = parent
                locus_primers.append(f)

            # try to find the same pcr product
            try:
                # if the position of the primers is equal, the sequences are 
                # also equal, so the termodynamic parameters should be the same as well
                # so we can use the Interval() '==' operator, which ignores gff attributes
                # the == operator raises error when the intervals contain each other
                def intervals_eq(a, b):
                    try:
                        return a == b
                    # when intervals contain each other, they're not equal
                    except NotImplementedError:
                        return False
            
                # find function with custom comparator
                def find(seq, item, _cmp):
                    for idx, i in enumerate(seq):
                        if _cmp(i, item): return idx
                    raise ValueError('item not found')

                idx = find(locus_primers, curr_primers['pcr-product'], intervals_eq)
                parent_pcr = locus_primers[idx]
            # new pcr product
            # add the product and pcr primers
            except ValueError:
                pcr_id += 1
                parent_pcr = curr_primers['pcr-product']
                parent_pcr.name = "%s-%d" % (locus_out_name, pcr_id)
                locus_primers.append(curr_primers['pcr-product'])
                
                # decorate the PCR product with reference genome transcript names 
                # and some statistics
                ref_names = set(e.attrs['Name'] for e in var_features if 'coords' in e.attrs)
                if len(ref_names): parent_pcr.attrs['ref_names'] = '|'.join(ref_names)
                parent_pcr.attrs['variants_in_locus'] = str(pass_variants)

                # add pcr primers
                add_feature(curr_primers['pcr-left'], "%sF" % parent_pcr.name, parent_pcr.name)
                add_feature(curr_primers['pcr-right'], "%sR" % parent_pcr.name, parent_pcr.name)

            # add genotyping primers
            if 'gt-left' in curr_primers:
                add_feature(curr_primers['gt-left'], "%s-%dF" % (parent_pcr.name, lvar_num + 1), parent_pcr.name)
            if 'gt-right' in curr_primers:
                add_feature(curr_primers['gt-right'], "%s-%dR" % (parent_pcr.name, lvar_num + 1), parent_pcr.name)
        
        # output all primers for current locus
        map(lambda f: sys.stdout.write(str(f)), locus_primers)

if __name__ == "__main__": main()
