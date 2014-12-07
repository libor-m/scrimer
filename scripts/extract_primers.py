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
        sys.exit('use: %s primer_gff [isPcr|table]' % sys.argv[0])

    # output blat isPcr format?
    out_type = 'fasta'
    if 'isPcr' in sys.argv:
        out_type = 'isPcr'
    if 'table' in sys.argv:
        out_type = 'table'

    pcr_pairs = dict()
    pcr_products = dict()
    # hack fix an error in newer pybedtools
    bt = pybedtools.BedTool(sys.argv[1])
    bt._isbam = False
    for f in bt:
        # fasta output
        if out_type == 'fasta':
            if 'SEQUENCE' in f.attrs:
                print ">", f.name
                print f.attrs['SEQUENCE']
        # isPcr output
        elif out_type == 'isPcr':
            if f.fields[2] != 'pcr-primer': continue
            
            pair = pcr_pairs.setdefault(f.attrs['Parent'], {})
            side = f.name[-1]
            pair[side] = f.attrs['SEQUENCE']
            if(len(pair) == 2):
                print f.name[:-1], pair['F'], pair['R']
                del pcr_pairs[f.attrs['Parent']]
        # descriptive table output
        # output only pcr pair and gt primer, if there is at least one valid
        # gt primer
        # the easiest way is to do this in two passes
        # - first fill a table with info on prc-products
        # - second output only the interesting stuff
        # this could be done in one pass using tabix index
        elif out_type == 'table':
            if f.fields[2] == 'rejected-var': continue

            if f.fields[2] == 'pcr-product':
                key = f.name
            else:
                key = f.attrs['Parent']

            pdata = pcr_products.setdefault(key, {})
            if f.fields[2] == 'pcr-product':
                pdata['product'] = f
            else:
                l = pdata.setdefault(f.fields[2], [])
                l.append(f)

    # output table
    # format:
    # variant-name pcr-f pcr-r gt-f gt-r product-size
    if out_type == 'table':
        fnames = ['ID', 'ref-chromosome', 'pcr-F', 'pcr-R', 'product-size', 'pcr-penalty', 'gt-tag',]
        # dp4 is: reference forward, reference reverse, alternate forward, alternate reverse
        att_names = ['SEQUENCE', 'TM', 'GC_PERCENT', 'PENALTY', 'VAR_mindps', 'VAR_dp4', 'VAR_FQ', 'VAR_MQ',]
        
        # print header
        print "\t".join(fnames + att_names)

        for id in sorted(pcr_products.iterkeys()):
            product = pcr_products[id]
            
            # missing gt primers?
            if 'gt-primer' not in product:
                continue
            
            # no valid gt primers
            if all('PROBLEMS' in f.attrs for f in product['gt-primer']):
                continue

            # get primers with orientations
            pcr = {f.name[-1]:f.attrs['SEQUENCE'] for f in product['pcr-primer']}
            
            # there can be several variants in single pcr product
            # duplicate the pcr product in the table for each possible genotyping primer
            gts = [(f.name[-2:], f) for f in product['gt-primer'] if 'PROBLEMS' not in f.attrs]
            
            for tag, f_gt in gts:
                fields = [id, 
                        product['product'].chrom,
                        pcr['F'], 
                        pcr['R'], 
                        product['product'].attrs['PRODUCT_SIZE'],
                        product['product'].attrs['PENALTY'],
                        tag, 
                        ]

                # add selected attributes
                fields.extend(f_gt.attrs[k] if k in f_gt.attrs else '' for k in att_names)

                print '\t'.join(fields)
            
if __name__ == "__main__": main()
