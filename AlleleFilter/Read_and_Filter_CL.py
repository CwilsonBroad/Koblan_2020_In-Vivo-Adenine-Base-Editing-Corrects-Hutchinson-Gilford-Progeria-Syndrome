#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 10:50:40 2021

@author: cwilson
"""

# ============================================================================
# Modules 
# ============================================================================ 
import sys

# ============================================================================
# Input Files 
# ============================================================================ 
Motif_File = 'Remove_List.txt'
Fastq_File = sys.argv[1]

# ============================================================================
# Read Motif File
# ============================================================================ 
Motif_List = []
with open('Remove_List.txt', 'r') as fh: 
    for line in fh.readlines(): 
        Motif_List.append(line[:-1])
        
# ============================================================================
# Read Fastq and Organize into Data
# ============================================================================
Fastq_List, Fastq_Entry = [], []
with open(Fastq_File, 'r') as fh:
    for ent, line in enumerate(fh.readlines()): 
        if ent % 4 == 0 and  ent != 0: 
            Fastq_List.append(Fastq_Entry)
            Fastq_Entry = []
        Fastq_Entry.append(line)
Fastq_List.append(Fastq_Entry)
     
# ============================================================================
# Organize Data into those with and without the motifs
# ============================================================================
Saved_List, Removed_List = [], []
for Fastq in Fastq_List: 
    Seq = Fastq[1][:-1]
    s = 1
    for Motif in Motif_List: 
        if Motif in Seq: 
            s = 0 
    if s == 1: 
        Saved_List.append(Fastq)
    else: 
        Removed_List.append(Fastq)

# ============================================================================
# Write Out Data
# ============================================================================
with open(Fastq_File[:-6] + '_Filtered.fastq', 'w') as out: 
    for Seq in Saved_List: 
        for i in Seq: 
            out.write(i)

with open(Fastq_File[:-6] + '_Discarded.fastq', 'w') as out: 
    for Seq in Removed_List: 
        for i in Seq: 
            out.write(i)

