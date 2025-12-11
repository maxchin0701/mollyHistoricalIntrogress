#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 15:16:09 2025

@author: maxchin
"""

#import module
from Bio import SeqIO
import argparse

#parse user inputs
parser = argparse.ArgumentParser(description="Splitting scaffolds into intervals of given length")
parser.add_argument("--inFasta", action="store", dest="inFile", type=str)
parser.add_argument("--scaff", action="store", dest="scaffs", type=str)
parser.add_argument("--dipProb", action="store", dest="dipProb", type=int)
parser.add_argument("--tripProb", action="store", dest="tripProb", type=int)
parser.add_argument("--outTSV", action="store", dest="outFile", type=str)
args=parser.parse_args()

#variables to store scaff names and length
seqN = []

#keep scaffolds above 15Mb
for record in SeqIO.parse(args.inFile, "fasta"):
    if record.name == args.scaff :
        seqN.append(record.name)
                
#write scaffolds to file       
out = open(args.outFile, 'w')
out.write("CONTIG_NAME\tPLOIDY_PRIOR_2\tPLOIDY_PRIOR_3\n")
for i in range(0,len(seqN)):
    out.write(seqN[i] + "\t" + str(args.dipProb) + "\t" + str(args.tripProb) + "\n")

out.close()


