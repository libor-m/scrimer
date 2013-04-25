#!/usr/bin/env python

"""Given pairs of 'almost identical' sequences create clusters of sequences.
From each cluster select the longest sequence. Output names of all other 
sequences (for example to remove them from the data afterwards).

Uses disjoint sets forest to store the clusters so it should scale to millions of sequences.

Input

- custom formated (``--format=general:name1,size1,start1,name2,size2,start2,strand2,identity,coverage``)
  output from ``lastz`` on standard input

Output 

- names of sequences that were selected as reduntant
"""

# This file is a part of Scrimer.
# See LICENSE.txt for details on licensing.
#    Copyright (C) 2012, 2013 Libor Morkovsky


import sys

class DisjointSet:
    """implementation by John Machin found at
    http://stackoverflow.com/questions/3067529/looking-for-a-set-union-find-algorithm
    """
    def __init__(self):
        self.leader = {} # maps a member to the group's leader
        self.group = {} # maps a group leader to the group (which is a set)

    def add(self, a, b):
        leadera = self.leader.get(a)
        leaderb = self.leader.get(b)
        if leadera is not None:
            if leaderb is not None:
                if leadera == leaderb: return # nothing to do
                groupa = self.group[leadera]
                groupb = self.group[leaderb]
                if len(groupa) < len(groupb):
                    a, leadera, groupa, b, leaderb, groupb = b, leaderb, groupb, a, leadera, groupa
                groupa |= groupb
                del self.group[leaderb]
                for k in groupb:
                    self.leader[k] = leadera
            else:
                self.group[leadera].add(b)
                self.leader[b] = leadera
        else:
            if leaderb is not None:
                self.group[leaderb].add(a)
                self.leader[a] = leaderb
            else:
                self.leader[a] = self.leader[b] = a
                self.group[a] = set([a, b])

def main():
    ds = DisjointSet()
    seqlens = {}
    
    # create clusters
    for line in sys.stdin:
        fields = line.split()
        
        # count on seqlen being always the same for same seq id
        def addlen(id, idlen):
            if id not in seqlens:
                seqlens[id] = int(idlen)
        
        addlen(fields[0], fields[1])
        addlen(fields[3], fields[4])
        
        ds.add(fields[0], fields[3])

    # find the longest sequence 
    # (for each item in each group find its length, store it to a sorted list)
    clusters = [sorted((seqlens[id], id) for id in g) for g in ds.group.itervalues()]

    # from pprint import pprint
    # pprint(clusters)

    # output all but the biggest sequence id
    for cluster in clusters:
        for _, id in cluster[:-1]:
            print id
    
    print >> sys.stderr, len(seqlens), "sequences in", len(clusters), "clusters"

if __name__ == "__main__": main() 