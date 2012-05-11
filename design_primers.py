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
# use the multiple records in single call to test different options
def find_pcr_primers(pseq, target):
    """call primer3 executable, return the results
    """
    def_params = {
    'PRIMER_THERMODYNAMIC_PARAMETERS_PATH':'/opt/primer3/bin/primer3_config/',
    'PRIMER_MAX_NS_ACCEPTED':'0',
    }
    
    p3 = primer3_connector.Primer3(**def_params)

    rec_params = {
    'SEQUENCE_ID':'design_primers_task',
    'SEQUENCE_TEMPLATE':pseq,
    'SEQUENCE_TARGET': "%d,1" % target,
    'PRIMER_OPT_SIZE':'19',
    'PRIMER_MIN_SIZE':'17',
    'PRIMER_MAX_SIZE':'25',
    'PRIMER_PRODUCT_SIZE_RANGE':'70-300',
    }
    
    return p3.call([rec_params])

#TODO: use some general score based on coverage, 
# exon support etc in the gff score field
def primer_to_gff(primer, tag, seq_start):
    pos, len = map(int, primer['position'].split(','))
    
    # transfer the calculated values to attributes
    # skip the fields used elsewhere in gff
    at = pybedtools.Attributes()
    for k, v in primer.iteritems():
        if k == 'PENALTY' or k == 'position': continue
        at[k] = v
        
    gflist = [var.CHROM, 'design-primers', tag, 
        seq_start + pos, seq_start + pos + len,
        primer['PENALTY'], '+', '.', str(at)]
        
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
def main():
    if len(sys.argv) < 2:
        sys.exit('use: %s genome_fasta_with_fai mrna_gff_with_tbi variants_vcf_with_tbi' % sys.argv[0])

    genome = pysam.Fastafile(sys.argv[1])
    annotations = pybedtools.BedTool(sys.argv[2])
    variants = vcf.Reader(filename=sys.argv[3])

    # use iterator filter, so the variants are streamed
    for var in itertools.ifilter(lambda v: not v.FILTER, variants):
        var_features = annotations.tabix_intervals(pybedtools.Interval(var.CHROM, var.POS, var.POS))
        var_exons = [f for f in var_features if f.fields[2] == 'exon']
        
        # get a minimal predicted exon
        # (trying to be conservative)
        # consider only predicted, not transferred exons (source in attrs)
        ex_start = max(f.start for f in var_exons if 'source' in f.attrs)
        ex_end = min(f.end for f in var_exons if 'source' in f.attrs)
        min_exon = pybedtools.Interval(var.CHROM, ex_start, ex_end)
        
        # check if the exon has a sane size (> 70)
        if min_exon.length < 70: continue

        # check if the variant is in position that permits PCR amplification
        # (more than 20 bases from both ends)
        # this implies a check if the variant is still inside of the minimal exon
        distances = (var.POS - min_exon.start, min_exon.end - var.POS)
        if min(distances) < 20: continue
        
        # get all variants for the minimal exon
        exon_variants = variants.fetch(var.CHROM, min_exon.start, min_exon.end)
        
        # get sequence for the minimal exon
        # and patch all the variants with Ns
        seq = genome.fetch(min_exon.chrom, min_exon.start, min_exon.end)
        pseq = patched_sequence(seq, min_exon.start, exon_variants)
        
        # check if there is at least 20 fixed bases on either side
        # of the target variant because of the genotyping primer
        gt_primers = get_flanking(pseq, var.POS, 25)
        if all(map(lambda s: 'N' in s, gt_primers)): continue
        
        # call primer3 to design PCR primers
        # mark only the variant of interest as a target
        #TODO: call primer3 only once with multiple records?
        primers = find_pcr_primers(pseq, var.POS - 1 - min_exon.start)
        
        # if something was found, output the first pair
        # (should be the best one)
        if len(primers[0]['PAIR']):
            print primer_to_gff(primers[0]['LEFT'][0], 'primer-pcr-left')
            print primer_to_gff(primers[0]['RIGHT'][0], 'primer-pcr-right')
        else:
            continue
        
        # if we've got PCR primers, try to check the genotyping ones

if __name__ == "__main__": main()


