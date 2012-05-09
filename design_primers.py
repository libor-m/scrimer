#!/usr/bin/env python

# given the
# - reference sequence (fasta/samtools fai index)
# - annotations (gff, has to contain mRNA and exon entries)
# - filtered variants (vcf, variants with PASS are used for the primer design)
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
# Author: Libor Morkovsky 2012

import sys

import pysam # FastaFile
import vcf
import pybedtools
from primer3_connector import primer3

record = primer3_connector.BoulderIO.parse(
"""SEQUENCE_ID=example
SEQUENCE_TEMPLATE=GTAGTCAGTAGACGATGACTACTGACGATGCAGACNACACACACACACACAGCACACAGGTATTAGTGGGCCATTCGATCCCGACCCAAATCGATAGCTACGATGACG
SEQUENCE_TARGET=37,21
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_OPT_SIZE=18
PRIMER_MIN_SIZE=15
PRIMER_MAX_SIZE=21
PRIMER_MAX_NS_ACCEPTED=3
PRIMER_PRODUCT_SIZE_RANGE=50-100
""")

record_no_res = primer3_connector.BoulderIO.parse(
"""SEQUENCE_ID=example
SEQUENCE_TEMPLATE=GTAGTCAGTAGACNATGACNACTGACGATGCAGACNACACACACACACACAGCACACAGGTATTAGTGGGCCATTCGATCCCGACCCAAATCGATAGCTACGATGACG
SEQUENCE_TARGET=37,21
PRIMER_TASK=pick_detection_primers
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=18
PRIMER_MIN_SIZE=15
PRIMER_MAX_SIZE=21
PRIMER_MAX_NS_ACCEPTED=1
PRIMER_PRODUCT_SIZE_RANGE=75-100
P3_FILE_FLAG=1
SEQUENCE_INTERNAL_EXCLUDED_REGION=37,21
""")

default_params = primer3_connector.BoulderIO.parse(
"""PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/opt/primer3/bin/primer3_config/
PRIMER_MAX_NS_ACCEPTED=0
PRIMER_EXPLAIN_FLAG=1
""")[0]

p3 = primer3_connector.Primer3(**default_params)

# test for single record
p3.call(record)
# test for multiple records
p3.call(record * 2)

# old egglib version
import egglib
import pysam
ff = pysam.Fastafile("33-virtual-genome/lx3.fasta")
p3 = egglib.wrappers.Primer3(ff.fetch('chr1', 0, 200), SEQUENCE_TARGET='30,20', PRIMER_MAX_NS_ACCEPTED=0, PRIMER_THERMODYNAMIC_PARAMETERS_PATH='/opt/primer3/bin/primer3_config/')
p3.find_primers()

PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/opt/primer3/bin/primer3_config/