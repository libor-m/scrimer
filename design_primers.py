#!/usr/bin/env python

# given the
# - reference sequence (fasta with samtools fai index)
# - annotations (gff3, has to contain exon entries)
# - filtered variants (vcf, primers are designed for variants with PASS)
# produces 
# - PCR and genotyping primers selected using primer3 (gff3)
# algorithm
# - there is only a few selected variants, so the least amount of work
#   will be to do the work only for variants
# - for each of the selected variants
#   - request exons
#   - apply the technical constraints 
#     (minimal primer length of 20 form the edge of an exon)
#   - patch exon sequence to mark positions of known variants
#   - find suitable genotyping primers
#   - design PCR primers to flank the (usable) genotyping primers
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
import primer3_connector

min_gt_primer_len = 20
max_gt_primer_len = 28

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
    'PRIMER_THERMODYNAMIC_PARAMETERS_PATH':'/opt/primer3/bin/primer3_config/',
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
    'PRIMER_THERMODYNAMIC_PARAMETERS_PATH':'/opt/primer3/bin/primer3_config/',
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
    'PRIMER_OPT_SIZE':'22',
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

def check_gt_primers(pseq, left_primer, right_primer_rc, **kwargs):
    """call primer3 executable, return the results
    """
    def_params = {
    'PRIMER_THERMODYNAMIC_PARAMETERS_PATH':'/opt/primer3/bin/primer3_config/',
    'PRIMER_MAX_NS_ACCEPTED':'0',
    'PRIMER_EXPLAIN_FLAG':'1', 
    'PRIMER_PICK_ANYWAY':'1',
    }
    
    def_params.update(kwargs)
    
    p3 = primer3_connector.Primer3(**def_params)

    rec_params = {
    'SEQUENCE_ID':'check_gt_primers_task',
    'SEQUENCE_TEMPLATE': pseq,
    'PRIMER_TASK':'check_primers',
    'PRIMER_MIN_TM':'50',
    'PRIMER_PRODUCT_SIZE_RANGE':'40-100',
    }
    
    if len(left_primer):
        rec_params['SEQUENCE_PRIMER'] = left_primer
    if len(right_primer_rc):
        rec_params['SEQUENCE_PRIMER_REVCOMP'] = right_primer_rc
    
    return p3.call([rec_params])

def primer_to_gff(name, primer, tag, seq_name, seq_start, strand, **kwargs):
    """Create a gff feature from 
    partially parsed primer3 results.
    """
    pos, len = map(int, primer['position'].split(','))
    
    # transfer the calculated values to attributes
    # skip the fields used elsewhere in gff
    at = pybedtools.Attributes()
    for k, v in primer.iteritems():
        if k == 'position': continue
        at[k] = v.replace(';', '%3B')
        
    at['ID'] = name
    
    # pass all optional params to attributes
    at.update(kwargs)
    
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

    # remove the pybedtools newline, that does not work well with print
    return str(pybedtools.create_interval_from_list(gflist)).strip()

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
    
def main():
    if len(sys.argv) < 2:
        sys.exit('use: %s genome_fasta_with_fai gff_with_tbi variants_vcf_with_tbi' % sys.argv[0])

    genome = pysam.Fastafile(sys.argv[1])
    annotations = pybedtools.BedTool(sys.argv[2])
    variants = vcf.Reader(filename=sys.argv[3])
    # open the variants once more, so the ifilter is not broken by .fetch
    var_fetcher = vcf.Reader(filename=sys.argv[3])
    
    primer_names = dict()

    # use iterator filter, so the variants are streamed
    for var in itertools.ifilter(lambda v: not v.FILTER, variants):
        var_features = annotations.tabix_intervals(pybedtools.Interval(var.CHROM, var.POS, var.POS))
        var_exons = [f for f in var_features if f.fields[2] == 'exon']
        
        def report_rejected_var(reason):
            """ Output a false feature holding the 
            reason for rejection of selected variant in 
            the design process.
            """
            at = pybedtools.Attributes()
            at['Note'] = reason
            at['color'] = "#bb0000"
            
            gflist = [var.CHROM, 'design-primers', 'rejected-var', 
                str(var.POS - max_gt_primer_len), str(var.POS + max_gt_primer_len),
                '0', '+', '.', str(at)]
                
            print str(pybedtools.create_interval_from_list(gflist)).strip()

        if len(var_exons) == 0: 
            report_rejected_var('no predicted exons')
            continue
        
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
