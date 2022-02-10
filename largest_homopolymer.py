#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#########################################################################
##                                                                     ##
##     largest_homopolymer.py                                          ##
##                                                                     ##
##     Gael A. Millot                                                  ##
##     Bioinformatics and Biostatistics Hub                            ##
##     Computational Biology Department                                ##
##     Institut Pasteur Paris                                          ##
##                                                                     ##
#########################################################################




################################ Aim


################################ End Aim


################################ Introduction


################################ End Introduction


################################ Acknowlegments


################################ End Acknowlegments


################################ Initialization


from os.path import exists # to check if file exists
import os, sys, glob, re, time # import sys for having sys.arg (see below), time for having time (to chrono the running), by default in the base distribution of python
import csv # allow to import cvs files (parse the lines correctly), by default in the base distribution of python
import pandas as pd # allow to export tsv files
from Bio import SeqIO # to chack that file is fasta


################################ End Initialization


################################ Parameters that need to be set by the user


################################ End Parameters that need to be set by the user


################################ Config import


input_path = sys.argv[1]  # 1st argument: pedigree file name, 1 ID par ligne, sys.argv takes arguments from the bash line command when running a .py script
output_path = open(sys.argv[2],"w")  # create (override) an output file using the 4th argument of the running command. See https://docs.python.org/3/library/functions.html?highlight=open#open


################################ End Config import

################################ Test


# input_path = "C:/Users/Gael/Documents/Git_projects/largest_homopolymer/dataset/integrases.fasta"
# output_path = "C:/Users/Gael/Documents/Git_projects/largest_homopolymer/res.tsv"


################################ End Test

################################ Recording of the initial parameters





################################ End Recording of the initial parameters


################################ Functions


def homopoly_detect(sequence):
    '''
    AIM
        Detect the largest homopolymer in a batch of DNA sequences
    WARNINGS
    ARGUMENTS
        sequence: a sequence formed of any kind of letters
    RETURN
        a data frame containing:
    REQUIRED PACKAGES
        None
    EXAMPLE
        homopoly_detect(sequence = seq)
    '''

    homopoly_max = -1
    nuc_max_list = []
    pos_max_list = []
    homopoly_max_list = []
    count = 0
    prev_nuc = None
    ini_pos = -1

    for pos, nuc in enumerate(sequence): # https://docs.python.org/3/library/functions.html#enumerate # warning, first nuc is position 0
        if nuc == prev_nuc:
            count += 1
        else:
            if homopoly_max < count & count >= 2:
                nuc_max = prev_nuc
                pos_max = ini_pos
                homopoly_max = count
                nuc_max_list = [nuc_max] # list is recreated
                pos_max_list = [pos_max]
                homopoly_max_list = [homopoly_max]
            elif homopoly_max == count & count >= 2:
                nuc_max_list.append(prev_nuc) # list is appended
                pos_max_list.append(ini_pos)
                homopoly_max_list.append(count)
            else:
                nothing = None

            count = 1
            ini_pos = pos + 1 # because pos starts at zero

        prev_nuc = nuc

    return nuc_max_list, [str(x) for x in pos_max_list], [str(x) for x in homopoly_max_list] #convert all into strings


################################ End Functions


################################ Pre-ignition checking


# reserved words
# end reserved words
# argument primary checking


file_exists = exists(input_path)
if not file_exists:
    print("\n\n================\n\nERROR: FILE DOES NOT EXISTS AT THIS PATH:\n")
    print(input_path)
    print("\n\n================\n\n")
    quit()


# end argument primary checking
# second round of checking and data preparation
# management of NA arguments
# end management of NA arguments
# management of NULL arguments
# end management of NULL arguments
# code that protects set.seed() in the global environment
# end code that protects set.seed() in the global environment
# warning initiation
# warn.count <- 0 # not required
# end warning initiation
# other checkings


with open(input_path, "r") as fasta:
    res = SeqIO.parse(fasta, "fasta") # False when `fasta` is empty, i.e. wasn't a FASTA file
    if not res:
        print("\n\n================\n\nERROR: FILE IS NOT A FASTA FILE:\n")
        print(input_path)
        print("\n\n================\n\n")
        quit()


# end other checkings
# reserved word checking
# end reserved word checking
# end second round of checking and data preparation
# package checking
# end package checking


################################ End pre-ignition checking


################################ Main code


################ Ignition

name = []
nucleotid = []
position = []
size = []

with open(input_path, 'r') as fasta: # open the fasta file, 'r' means read, advantage of with command: cleanly close a file when the block ends (keep memory and avoid using close at the end). Check if and avoid error if 'w' writing option is selected) ?
    for line in fasta:
        if line.startswith(">"):
            name.append(line.strip()) # strip() added to remove the \n added at each end on "line" when appended in name
        else:
            nuc_max_list, pos_max_list, homopoly_max_list = homopoly_detect(line)
            if len(nuc_max_list) > 1:
                nuc_max_list = ";".join(nuc_max_list)
                pos_max_list = ";".join(pos_max_list)
                homopoly_max_list = ";".join(homopoly_max_list)
            else:
                nuc_max_list = "".join(nuc_max_list)
                pos_max_list = "".join(pos_max_list)
                homopoly_max_list = "".join(homopoly_max_list)
            nucleotid.append(nuc_max_list)
            position.append(pos_max_list)
            size.append(homopoly_max_list)


df = pd.DataFrame(list(zip(name, nucleotid, position, size)), columns =['name', 'nucleotide', 'position', 'size'])
df.to_csv(path_or_buf = output_path, sep = "\t", index = False)


################ End ignition


################ Graphical parameter initialization


################ End graphical parameter initialization


################ Data import


################ end Data import


############ modifications of imported tables


############ end modifications of imported tables


############ plotting


############ end plotting


################ Pdf window closing


################ end Pdf window closing


################ Seeding inactivation


################ end Seeding inactivation


################ Environment saving


################ end Environment saving


################ Warning messages


################ end Warning messages


################ Parameter printing


################ end Parameter printing


################################ End Main code







