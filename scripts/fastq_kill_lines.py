#!/usr/bin/env python

"""
Input

- FASTQ file 
- list of indices (0 based)
 
Output

- FASTQ file without the given sequences 

Author: Libor Morkovsky, 2012
"""

# This file is a part of Scrimer.
# See LICENSE.txt for details on licensing.
#    Copyright (C) 2012, 2013 Libor Morkovsky

import sys

def main():
    # read the kill list
    fkill = open(sys.argv[1])
    killlines = fkill.readlines()
    killseq = frozenset([int(line) for line in killlines])
    fkill.close()

    # go through the input, avoid outputting the 4 line blocks 
    # with given indices
    fseq = open(sys.argv[2])
    linenum = 0
    for line in fseq:
        seqnum = linenum / 4
        linenum += 1  
        if seqnum in killseq:
            continue
        print line.strip()

    fseq.close()

if __name__ == "__main__": main()
