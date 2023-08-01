#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#########################################################################
##                                                                     ##
##     homopolymer.py                                                  ##
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


import sys # import sys for having sys.arg (see below)
import pandas as pd # allow to export tsv files
import numpy # for numpy.multiply()
import regex # for sub() function
import random # for shuffle() function


################################ End Initialization


################################ Parameters that need to be set by the user


################################ End Parameters that need to be set by the user


################################ Config import


seq = sys.argv[1]  # 1st argument: pedigree file name, 1 ID par ligne, sys.argv takes arguments from the bash line command when running a .py script
name = sys.argv[2]
min_length = sys.argv[3]
output_path = sys.argv[4]


################################ End Config import

################################ Test


# seq = "ATGAAATCCCAGTTTTTGTTAAGTGTTCGCGAATTTATGCAAACTCGATACTATGCAAAAAAAACCATAGAAGCTTACCTTCATTGGATCACTCGTTACATCCATTTTCATAATAAAAAGCACCCTAGCTTAATGGGAGATAAAGAGGTCGAAGAATTTTTAACCTACTTAGCCGTGCAAGGTAAAGTGGCCACAAAGACTCAATCACTAGCCTTGAACTCACTCAGTTTTCTATACAAAGAAATTCTAAAAACACCCCTTT"
# name = "lcl|NC_002506.1_cds_NP_232687.1_1"
# min_length = 2
# output_path = "C:/Users/Gael/Documents/Git_projects/homopolymer/res.tsv"



################################ End Test

################################ Recording of the initial parameters




################################ End Recording of the initial parameters


################################ Functions



def homopoly_detect(sequence, mini):
    '''
    AIM
        Detect the largest homopolymer in a batch of DNA sequences
    WARNINGS
    ARGUMENTS
        sequence: a sequence formed of any kind of letters
        mini: the mini length of homopolymer
    RETURN
        a list containing:
            size_max_list: list of the length of the longest homopol (list because of the case of equality, but here, all the numbers are the saze because same max size. Differences will be for the related nuc_max_list and pos_max_list)
            nuc_max_list: list of the nucleotid (A,T,C,G) that is the longest homopol (list because of the case of equality)
            pos_max_list: list of the position (using the first nuc of the homopol) that is the longest homopol (list because of the case of equality)
            homopoly_nb: number of homopolymers
            homo_distrib_list: list of all the number of homopol, the first position in the list being the homopol of size 1 and the last, the homopol of size of the sequence
    REQUIRED PACKAGES
        None
    EXAMPLE
        homopoly_detect(sequence = seq, mini = 2)
    DEBUGGING
        sequence = "ATT"
    '''
    
    length = len(sequence.strip())

    size_max = -1
    nuc_max_list = [] # list is recreated
    pos_max_list = []
    size_max_list = []
    count = 0
    prev_nuc = None
    homopoly_nb = 0 # nb of homopolymers in the seq
    # Do not fill data frames: first list then convert https://stackoverflow.com/questions/13784192/creating-an-empty-pandas-dataframe-then-filling-it
    homo_distrib_list = [0] * length # list of zero that will contain the count of all the homopolymers, the first position in the list being the homopol of size 1 and the last, the homopol of size of the sequence

    for pos, nuc in enumerate(sequence): # https://docs.python.org/3/library/functions.html#enumerate # warning, first nuc is position 0 bt I could have written enumerate(sequence, start = 1) # split the sequence in single nuc and provide the pos of each of these nuc (0, 1, 2, etc.). Thus, pos will be 0, 1, 2, etc. and nuc will be the associated nuc of the seq (ex: A,T,G, etc.)
        if length == 1: # input seq of 1 nucleotide
            if mini == 1: # else -> the output lists remain empty
                nuc_max_list = [nuc]
                pos_max_list = [pos + 1]
                size_max_list = [1]
                homopoly_nb += 1
                homo_distrib_list[0] += 1
            
        else:
            if pos + 1 == 1: # first nucleotid
                count = 1
                nuc_max = nuc
                pos_max = pos + 1
                size_max = count
            
            else: # between second and last nuc included
                if nuc == prev_nuc:
                    count += 1
                    
                if nuc != prev_nuc: # warning: even if at pos, we are taking the info of prev_nuc here
                    homo_distrib_list[count - 1] += 1
                    if count >= mini:
                        homopoly_nb += 1
                        if size_max < count:
                            nuc_max = prev_nuc
                            pos_max = pos + 1 - count
                            size_max = count
                            nuc_max_list = [nuc_max] # list is recreated
                            pos_max_list = [pos_max]
                            size_max_list = [size_max]

                        elif size_max == count:
                           nuc_max_list.append(prev_nuc) # list is appended
                           pos_max_list.append(pos + 1 - count)
                           size_max_list.append(size_max)
                        
                        else:
                            nothing = None
            
                    count = 1 # reset counting
                    
                if pos + 1 == length: # last nucleotid
                    homo_distrib_list[count - 1] += 1
                    if count >= mini:
                        homopoly_nb += 1
                        if size_max < count:
                            nuc_max = nuc
                            pos_max = pos + 2 - count
                            size_max = count
                            nuc_max_list = [nuc_max] # list is recreated
                            pos_max_list = [pos_max]
                            size_max_list = [size_max]

                        elif size_max == count:
                           nuc_max_list.append(nuc) # list is appended
                           pos_max_list.append(pos + 2 - count)
                           size_max_list.append(size_max)
                        
                        else:
                            nothing = None
        prev_nuc = nuc

    return size_max_list, nuc_max_list, pos_max_list, homopoly_nb, homo_distrib_list # warning: homopoly_nb is not a list



################################ End Functions


################################ Pre-ignition checking


# reserved words
# end reserved words
# argument primary checking
min_length=int(min_length)
# end argument primary checking
# second round of checking and data preparation
# management of NA arguments
# end management of NA arguments
# management of NULL arguments
# end management of NULL arguments
# code that protects set.seed() in the global environment
random.seed(1)
# end code that protects set.seed() in the global environment
# warning initiation
# warn.count <- 0 # not required
# end warning initiation
# other checkings
# end other checkings
# reserved word checking
# end reserved word checking
# end second round of checking and data preparation
# package checking
# end package checking


################################ End pre-ignition checking


################################ Main code


################ Ignition



# observed homopolymers
length = len(seq.strip()) # length of the string (here the sequence)
size_max_list, nuc_max_list, pos_max_list, homopoly_nb, homo_distrib_list = homopoly_detect(sequence = seq, mini = min_length)
tempo_verif = numpy.multiply(homo_distrib_list, list(range(1, len(homo_distrib_list) + 1))).sum() #  + 1 to deal with the fact that range exclude the last number
if length != tempo_verif:
    sys.exit("Error in homopolymer.py: the observed distribution list homo_distrib_list should have a total sum corresponding to the length of the sequence")

homopoly_rel_pos_list = [(x - 1) / (length) for x in pos_max_list]
pos_max_list = [str(x) for x in pos_max_list] # to convert list of integers into strings
size_max_list = [str(x) for x in size_max_list] # to convert list of integers into strings
homopoly_rel_pos_list = [str(x) for x in homopoly_rel_pos_list] # to convert list of integers into strings
if len(nuc_max_list) > 1:
    nuc_max_list = ";".join(nuc_max_list) # combine (collapse) with ; as separator
    pos_max_list = ";".join(pos_max_list)
    size_max_list = ";".join(size_max_list)
    homopoly_rel_pos_list = ";".join(homopoly_rel_pos_list)
else:
    nuc_max_list = "".join(nuc_max_list)
    pos_max_list = "".join(pos_max_list)
    size_max_list = "".join(size_max_list)
    homopoly_rel_pos_list = "".join(homopoly_rel_pos_list)
# end observed homopolymers


# random homopolymers
n = 10000 # number of random samplings
sum_rd_homo_distrib_list = [0] * length
for i0 in list(range(1, n)):
    tempo_seq=''.join(random.sample(seq, len(seq))) # shuffle the single nuc of seq and return a list, which is then joined
    x1, x2, x3, x4, tempo_homo_distrib_list = homopoly_detect(sequence = tempo_seq, mini = 1)  # Warning: here min = 1 and not min_length because we only want the distrib of the homopolymers
    sum_rd_homo_distrib_list = [x + y for x, y in zip(sum_rd_homo_distrib_list, tempo_homo_distrib_list)]

mean_rd_homo_distrib_list = [x / n for x in sum_rd_homo_distrib_list]
mean_rd_homo_distrib_list = [str(x) for x in mean_rd_homo_distrib_list] # to convert list of integers into strings
# end random homopolymers


# data frame creation for export
tempo_num = numpy.multiply(homo_distrib_list[(min_length - 1):(len(homo_distrib_list) + 1)], list(range(min_length, len(homo_distrib_list) + 1))).sum() #  [1:3] starts at 0 and exclude the last number take second and third position of the list, excluding the fourth (3)
tempo_denom = numpy.array(homo_distrib_list[(min_length - 1):(len(homo_distrib_list) + 1)]).sum()
mean_size = tempo_num / tempo_denom
homo_distrib_list = [str(x) for x in homo_distrib_list] # to convert list of integers into strings
homo_distrib_list = ";".join(homo_distrib_list)
mean_rd_homo_distrib_list = ";".join(mean_rd_homo_distrib_list)
df = pd.DataFrame(columns =['name', 'seq_length', 'max_size', 'nucleotide', 'starting_position', 'relative_position', "nb", 'mean_size', "homopol_freq_according_to_size", str(n) + "_mean_rd_freq_according_to_size"])
df.loc[0] = [name, length, size_max_list, nuc_max_list, pos_max_list, homopoly_rel_pos_list, homopoly_nb, mean_size, homo_distrib_list, mean_rd_homo_distrib_list]
# name: name of the sequence
# seq_length: length of the sequence
# nucleotid: nuc of the homopolymer
# starting_position: first base position of the homopolymer
# relative_position: relative position of the homopolymer in the sequence using  y = (x - 1) / (size - length) with x the first position of the homopolymer, size the sequence length and length the homopolymer length (0 <= y <= 1)
# size: length of the homopol
# nb: homopolymer nb in seq (considering polymers starting at length 1)
# mean_size: average size of homopolymers (total number of homopolymer sof size between min_length and seq_length. It corresponds to homopoly_nb / (numpy.array(homo_distrib_list[min_length:(len(homo_distrib_list) + 1)]).sum())
# homopol_freq_according_to_size: 
df.to_csv(path_or_buf = output_path, sep = "\t", index = False, header = False)
# end data frame creation for export

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
random.seed()

################ end Seeding inactivation


################ Environment saving


################ end Environment saving


################ Warning messages


################ end Warning messages


################ Parameter printing


################ end Parameter printing


################################ End Main code







