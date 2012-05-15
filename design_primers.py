#!/usr/bin/env python

# given the
# - reference sequence (fasta with samtools fai index)
# - annotations (gff, has to contain mRNA and exon entries)
# - filtered variants (vcf, primers are designed for variants with PASS)
# produces 
# - PCR primers selected using primer3 
#   and genotyping primers checked by primer3 (gff)
# algorithm
# - there is only a few selected variants, so the least amount of work
#   will be to do the work only for variants
# - for each of the selected variants
#   - request exons
#   - apply the technical constraints 
#     (minimal primer length of 20 form the edge of an exon)
#   - patch exon sequence to mark positions of known variants
#   - design PCR primers
#   - check genotyping primers
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

#TODO: add ability to pick worse PCR primers in exchange
# for having a better genotyping primer (or having it at all)
# i.e. if the PCR primer on one side is very close to the 
# variant site, but there is another variable site inside 
# the only remaining possible genotyping primer
# - use the multiple records in single call to test different options?
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
    'SEQUENCE_TARGET': "%d,1" % target,
    'PRIMER_OPT_SIZE':'19',
    'PRIMER_MIN_SIZE':'17',
    'PRIMER_MAX_SIZE':'25',
    'PRIMER_PRODUCT_SIZE_RANGE':'70-300',
    }
    
    return p3.call([rec_params])

# not used, documentation of primer3 is not clear enough
# about how to use this feature..
def find_gt_primers(pseq, target):
    """call primer3 executable, return the results
    """
    def_params = {
    'PRIMER_THERMODYNAMIC_PARAMETERS_PATH':'/opt/primer3/bin/primer3_config/',
    'PRIMER_MAX_NS_ACCEPTED':'0',
    }
    
    p3 = primer3_connector.Primer3(**def_params)

    rec_params = {
    'SEQUENCE_ID':'pick_gt_primers_task',
    'SEQUENCE_TEMPLATE':pseq,
    'SEQUENCE_TARGET': "%d,1" % target,
    'PRIMER_TASK':'pick_discriminative_primers',
    'SEQUENCE_INCLUDED_REGION':  "%d,1" % target,
    'PRIMER_OPT_SIZE':'22',
    'PRIMER_MIN_SIZE':'20',
    'PRIMER_MAX_SIZE':'25',
    'PRIMER_PRODUCT_SIZE_RANGE':'1-100',
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

#TODO: use some general score based on coverage, 
# exon support etc in the gff score field
def primer_to_gff(name, primer, tag, seq_name, seq_start, strand):
    pos, len = map(int, primer['position'].split(','))
    
    # transfer the calculated values to attributes
    # skip the fields used elsewhere in gff
    at = pybedtools.Attributes()
    for k, v in primer.iteritems():
        if k == 'PENALTY' or k == 'position': continue
        at[k] = v
        
    at['Name'] = name
    
    gflist = [seq_name, 'design-primers', tag, 
        str(seq_start + pos + 1), str(seq_start + pos + len),
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
    
def main():
    if len(sys.argv) < 2:
        sys.exit('use: %s genome_fasta_with_fai mrna_gff_with_tbi variants_vcf_with_tbi' % sys.argv[0])

    genome = pysam.Fastafile(sys.argv[1])
    annotations = pybedtools.BedTool(sys.argv[2])
    variants = vcf.Reader(filename=sys.argv[3])
    # open the variants once more, so the ifilter is not broken by .fetch
    var_fetcher = vcf.Reader(filename=sys.argv[3])

    # use iterator filter, so the variants are streamed
    for var in itertools.ifilter(lambda v: not v.FILTER, variants):
        var_features = annotations.tabix_intervals(pybedtools.Interval(var.CHROM, var.POS, var.POS))
        var_exons = [f for f in var_features if f.fields[2] == 'exon']
        
        def report_rejected_var(reason):
            print '#', var, 'was rejected:', reason

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
            report_rejected_var('variant too close to minimal exon boundary')
            continue
        
        # get all variants for the minimal exon
        exon_variants = var_fetcher.fetch(var.CHROM, min_exon.start, min_exon.end)
        
        # get sequence for the minimal exon
        # and patch all the variants with Ns
        seq = genome.fetch(min_exon.chrom, min_exon.start, min_exon.end)
        pseq = patched_sequence(seq, min_exon.start, exon_variants)
        
        # check if there is at least 20 fixed bases on either side
        # of the target variant because of the genotyping primer
        gt_primer_seqs = get_flanking(pseq, var.POS - 1 - min_exon.start, 25)
        # attach 'N' to end of each primer, so we get the full len
        # instead of -1, when there is no N in the original primer
        max_gt_lens = [(s + 'N').find('N') for s in  [gt_primer_seqs[0][::-1], gt_primer_seqs[1]]]
        if all(x < 20 for x in max_gt_lens): 
            report_rejected_var('no possible genotyping primer longer than 20 bp')
            continue
        
        # call primer3 to design PCR primers
        # mark only a single base - the variant of interest - as a target
        #TODO: call primer3 only once with multiple records?
        primers = find_pcr_primers(pseq, var.POS - 1 - min_exon.start)
        
        # generate the name for the primer ensemble
        # use the name from reference genome only if it is unique
        # otherwise use the transcript name
        ref_names = set(e.attrs['Name'] for e in var_exons if 'coords' in e.attrs)
        if len(ref_names) == 1:
            primer_name = ref_names.pop()
        else:
            names = [e.attrs['Name'] for e in ex_predicted if 'Name' in e.attrs]
            primer_name = names[0] if len(names) else 'NAME-UNKNOWN'
        
        # if something was found, output the first pair
        # (should be the best one)
        if len(primers[0]['PAIR']):
            print primer_to_gff('LU-PCR-%s-F' % primer_name, primers[0]['LEFT'][0], 'primer-pcr', var.CHROM, min_exon.start, '+')
            print primer_to_gff('LU-PCR-%s-R' % primer_name, primers[0]['RIGHT'][0], 'primer-pcr', var.CHROM, min_exon.start, '-')
        else:
            report_rejected_var('no suitable PCR primers found')
            continue
        
        #TODO: if we've got PCR primers, try to check the genotyping ones
        # BatchPrimer3 is checking 
        # - melting temperature - GC content - repeats - Ns
        # - self complementarity (by Smith-Waterman in perl)
        # primer3 can be used to check GC, Tm, hairpins, self complementarity by thermodynamic approach
        # that seems superior to perl-based approach
        # pass sequences cut at the first N in the 25 bp input
        gt_primers = check_gt_primers(pseq, gt_primer_seqs[0][-max_gt_lens[0]:], reverse_complement(gt_primer_seqs[1][:max_gt_lens[1]]))
        
        # this should always report some results, so output those to gff
        print primer_to_gff('LU-GT-%s-F' % primer_name, gt_primers[0]['LEFT'][0], 'primer-gt', var.CHROM, min_exon.start, '+')
        print primer_to_gff('LU-GT-%s-R' % primer_name, gt_primers[0]['LEFT'][0], 'primer-gt', var.CHROM, min_exon.start, '-')
        
        #TODO: here should be some kind of scoring for the primer ensemble
        # - number of other segregating polymorphisms in the same exon/?mRNA 
        # - variant/primer coverage summary (would require .bam or some wiggle 
        #   file)

if __name__ == "__main__": main()
