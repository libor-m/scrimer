#!/usr/bin/env python

# given the
# - reference sequence (fasta with samtools fai index)
# - annotations (gff, has to contain mRNA and exon entries)
# - filtered variants (vcf, primers are designed for variants with PASS)
# produces 
# - PCR primers selected using primer3 
#   and genotyping primers checked by primer3 (gff)
# algorithm
# - go through the gff
# - for each mRNA get the variants and exons
#   go on if there are no variants or exons
# - for each unfiltered variant apply the technical constraints 
#   (minimal primer length of 20 form the edge of an exon)
# - construct a PCR candidate sequence from the exon containing the current 
#   variant of interest by substituting all variable sites with Ns
# algorithm2
# - there is only a few selected variants, so the least amount of work
#   will be to do the work only for variants
# - for each of the selected variants
#   - request exons
#   - apply the technical constraints 
#     (minimal primer length of 20 form the edge of an exon)
#   - 
# Author: Libor Morkovsky 2012
#

import sys
import itertools

import pysam # FastaFile
import vcf
import pybedtools
from primer3_connector import primer3


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
        #FIXME:unfinished
        
    for gffeature in annotations:
        # we're interested only in mRNA entries in the virtual genome
        # (those represent the original assembled contigs)
        if gffeature.fields[2] != 'mRNA': continue
        
        # get all exons for the mRNA
        contig_features = annotations.tabix_intervals(gffeature)
        exons = [f for f in contig_features if f.fields[2] == 'exon']
        
        
        

if __name__ == "__main__": main()


