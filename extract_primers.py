#!/usr/bin/env python

"""
Input

- gff file containing primers
- optional output format

Output

- fa/isPcr file with the primer sequences

Algorithm

- for each line in gff

 - **fa output**:
   if current line has a ``SEQUENCE`` attribute, output it
 - **ispcr output**:
   only pcr primers are of interest, and in equivalent pairs
   on encountering gt-xx put it to dict keyed by ``Parent``
   if there are both entries output and remove from dict

Author: Libor Morkovsky, 2012
"""

# This file is a part of Scrimer.
# See LICENSE.txt for details on licensing.
#    Copyright (C) 2012, 2013 Libor Morkovsky

import sys
import os
import itertools
import string

import pybedtools

    
def main():
    if len(sys.argv) < 2:
        sys.exit('use: %s primer_gff [isPcr]' % sys.argv[0])

    # output blat isPcr format?
    b_isPcr = 'isPcr' in sys.argv
    pcr_pairs = dict()
    for f in pybedtools.BedTool(sys.argv[1]):
        if not b_isPcr:
            if 'SEQUENCE' in f.attrs:
                print ">", f.name
                print f.attrs['SEQUENCE']
        else:
            if f.fields[2] != 'pcr-primer': continue
            
            pair = pcr_pairs.setdefault(f.attrs['Parent'], {})
            side = f.name[-1]
            pair[side] = f.attrs['SEQUENCE']
            if(len(pair) == 2):
                print f.name[:-1], pair['F'], pair['R']
                del pcr_pairs[f.attrs['Parent']]

if __name__ == "__main__": main()
