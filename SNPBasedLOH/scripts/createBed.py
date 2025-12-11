#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 13:12:44 2024

@author: maxchin
"""

#import module
from Bio import SeqIO
import argparse

#parse user inputs
parser = argparse.ArgumentParser(description="Splitting scaffolds into intervals of given length")
parser.add_argument("--inFasta", action="store", dest="inFile", type=str)
parser.add_argument("--scaff", action="store", dest="scaff", type=str)
parser.add_argument("--intSize", action="store", dest="intSize", type=int)
parser.add_argument("--outBed", action="store", dest="outFile", type=str)
args=parser.parse_args()

#variables to store scaff names and length
seqN = []
seqL = []

#keep specific scaffold
for record in SeqIO.parse(args.inFile, "fasta"):
    if record.name == args.scaff :
        seqL.append(len(record.seq))
        seqN.append(record.name)

#variables to store scaff names, start, and end
scaffL = []
startL = []
stopL = []
        
#loop through scaffolds
for i in range(0,len(seqN)) :
    
    #define name, start, and end for current scaffold
    scaffCur = seqN[i]
    startCur = 0
    stopCur = 0
    
    #while loop to iterate until we exceed scaffold length
    while stopCur < (seqL[i] - 1) :
        #append scaffold and current start
        scaffL.append(scaffCur)
        startL.append(startCur)
        
        #get current end
        stopCur = startCur + args.intSize
        
        #check if current stop is larger than or equal to scaffold size
        if stopCur >= seqL[i]:
            stopCur = seqL[i] - 1
        
        #append stopCur
        stopL.append(stopCur)
        
        #set startCur to new value
        startCur = startCur + args.intSize

#write file
print("Writing .bed file")
out = open(args.outFile, 'w')

for i in range(0,len(scaffL)):
    out.write(scaffL[i] + "\t" + str(startL[i]) + "\t" + str(stopL[i]) + "\n")

out.close()

print("Done")

