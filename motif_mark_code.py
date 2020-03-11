#!/bin/python

# Author: Wyatt Eng
# Date 2/17/20
# Goal: Parse a FASTA file that contains a gene (with introns and exons) and identify biological motifs. Motifs short sequences
#       of DNA that is indicative of structural elements of proteins. The input of the software should be a FASTA file and
#       a file containing the known motifs. The output should be a graphical version that identifies every motif location.

import cairo
import math
import argparse
import re
import os
from Bio import SeqIO
from matplotlib import cm
import random as random

# Set input variables
def get_args():
    parser = argparse.ArgumentParser(
        description="Identify motifs from gene files and produce a figure showing each motif's location across exons and introns. \
        Uppercase nucleotides should signify exons and Lowercase nucleotides should signify introns.")
    parser.add_argument("-f", "--file_name", help="path to a fasta file containing gene sequences", type=str)
    parser.add_argument("-m", "--motif_file_name", help="path to a file containing a list of the known motifs", type=str)
    parser.add_argument("-o", "--out_file_png_name", help="Desired outfile path + name for the output png file", type=str)
    parser.add_argument("-t", "--molecule_type", help="Molecule type of fasta (gene) file: RNA or DNA", type=str)
    return parser.parse_args()

args = get_args()
# Setting user inputted variables
input_file = args.file_name
motif_file = args.motif_file_name
out_file = args.out_file_png_name
fasta_convert = args.molecule_type

##### !!!!!DELETE LATER!!!!! #####
# motif_file = "Fig_1_motifs.txt"
# input_file = "Figure_1.fasta"
# out_file = "test.png"
# fasta_convert = "dna"

fasta_convert = fasta_convert.upper()

##### Global Variables ######
# Dictionary to store every possible motif
motifDict = {}
# Dictionary to associate a number to a motif
motif_numb = {}
# List storing the motif lengths
length_motif = []
# List storing the differences in lengths of the motifs
uniq_diff_motif_len = []
# Dictionary for amiguous nucleotides
ambig_dict = {'U':['T'], 'W':['A', 'T'], 'S':['C','G'], 'M':['A','C'], 'K':['G','T'], 'R':['A','G'], 'Y':['C','T'], 'B':['C','G','T'], 'D':['A','G','T'], 'H':['A','C','T'], 'V':['A','C','G'], 'N':['A','C','G','T']}

class cMotif:
  def __init__(self, mStart_pos, mStop_pos, mNumber):
      self.mStart_pos = mStart_pos
      self.mStop_pos = mStop_pos
      self.mColor = mNumber

class cGene:
  def __init__(self, mLength, mExonStart, mExonStop, mHeaderName, mMotifList, mGeneNum, mMotifDict):
      self.mLength = mLength
      self.mExonStart = mExonStart
      self.mExonStop = mExonStop
      self.mHeaderName = mHeaderName
      self.mMotifList = mMotifList
      self.mGeneNum = mGeneNum
      self.mMotifDict = mMotifDict

def f_IterateMotif(original, sequence, index):
    '''Function to iterate through every possible motif accounting for U and Y in the motif sequence. \
    Dictionary key: different possibilities of motifs / value: given motif \
    U (purines) = T or U / Y (pyrimidines) = C or T or U \
    Recursive Function! Shoutout to Clio Skevington) \
    Returns: Nothing, stores all possible motifs to dictionary'''
    if index == len(original):
        motifDict[sequence] = original
    elif sequence[index] in ambig_dict.keys():
        current_key = sequence[index]
        for nuc in ambig_dict[current_key]:
            seq_list = list(sequence)
            seq_list[index] = nuc
            sequence = "".join(seq_list)
            f_IterateMotif(original, sequence, index + 1)
    else:
        f_IterateMotif(original, sequence, index + 1)

def f_SaveMotifs(motif_file_name):
    '''Function to take in the motif file name and save known motifs in a dictionary. \
     Dictionary key: different possibilities of motifs / value: given motif \
     U (purines) = A or G / Y (pyrimidines) = C or T\
     Returns: Nothing'''
    index_val = 1
    with open(motif_file_name) as mf:
        for line in mf:
            line = line.strip()
            if line != "":
                line = line.upper()
                length_motif.append(len(line))
                f_IterateMotif(line, line, 0)
                motif_numb.setdefault(line, index_val)
                index_val = index_val + 1
    length_motif.sort(reverse=True)
    shortest_motif = length_motif[-1]
    diff_motif_len = []
    for i in range(0, (len(length_motif)-1)):
        diff_motif_len.append(abs(shortest_motif-length_motif[i]))
    diff_motif_len.sort()
    for x in diff_motif_len:
        # check if exists in unique_list or not
        if x not in uniq_diff_motif_len:
            uniq_diff_motif_len.append(x)


def multi2linefasta(file):
    '''Function converts multiline sequence fasta files to single line sequence fasta files \
    Returns: Nothing'''
    mfasta = "singleline.fasta"
    infile = open(file,'r')
    with open(mfasta, 'w') as ofile:
        for record in SeqIO.parse(infile, "fasta"):
            sequence = str(record.seq)
            if fasta_convert == "RNA":
                sequence = sequence.replace("u","t")
                sequence = sequence.replace("U", "T")
            ofile.write('>'+record.description+'\n'+sequence+'\n')
    infile.close()

def f_ParseSequence(sequence):
    '''Function takes in a sequence and finds the possible motifs in the sequence \
    Returns: Dictionary of motif (key: known motifs, value: number of appearances), list of motifs'''
    temp_dict = {}
    temp_motif_list = []
    with open(motif_file, "r") as mf:
        for line in mf:
            line = line.strip()
            line = line.upper()
            temp_dict.setdefault(line, 0)

    counter = 0
    for i in range(length_motif[-1], len(sequence)):
        sub_string = sequence[counter:i]
        sub_string = sub_string.upper()
        if sub_string in motifDict and len(sub_string) in length_motif:
            motif = motifDict.get(sub_string)
            numb = motif_numb.get(motif)
            temp_motif_list.append(cMotif(counter + 1, i, numb))
            motif_count = temp_dict[motif]
            temp_dict[motif] = motif_count + 1

        for val in uniq_diff_motif_len:
            sub_string = sequence[counter:i+val]
            sub_string = sub_string.upper()
            if sub_string in motifDict and len(sub_string) in length_motif:
                motif = motifDict.get(sub_string)
                numb = motif_numb.get(motif)
                temp_motif_list.append(cMotif(counter + 1, i, numb))
                motif_count = temp_dict[motif]
                temp_dict[motif] = motif_count + 1

        counter = counter + 1

    return temp_dict, temp_motif_list

def f_draw(Complete_Gene_List, lonest_gene):
    color_dict = {}
    color_count = 1
    numb_motifs = len(motif_numb)
    color_step = 255 / numb_motifs
    for key in motif_numb.keys():
        color_val = color_count * color_step
        tuple = cm.jet(int(color_val))
        color_dict.setdefault(key, tuple[0:3])
        color_count = color_count + 1

    temp_dict = {} # key = color number / value = motif
    for key in motif_numb.keys():
        temp_dict.setdefault(motif_numb[key], key)

    # Calculate Width
    longest_motif_font_length = 12 * length_motif[0]
    left_border = longest_motif_font_length + 37
    num_genes = len(Complete_Gene_List)
    height = num_genes * 150
    width = left_border + 30 + lonest_gene

    # Initialize PyCairo Surface
    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
    ct = cairo.Context(surface)

    # Draw a white background
    ct.set_source_rgb(1, 1, 1)
    ct.rectangle(0, 0, width, height)  # X1, Y1, Width, Height
    ct.fill()

    # Draw Legend
    motif_count = 0
    for key in color_dict.keys():
        # Draw Motif Color Boxes
        motif_height = motif_count * 15 + 20
        tup_color = color_dict[key]
        ct.set_source_rgb(tup_color[0], tup_color[1], tup_color[2])
        ct.set_line_width(3)
        ct.rectangle(10, motif_height, 15, 10)
        ct.fill()

        # Draw Motif Names
        ct.set_source_rgb(0.1, 0.1, 0.1)
        ct.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        ct.set_font_size(12)
        ct.move_to(28, motif_height + 10)
        ct.show_text(key)
        motif_count = motif_count + 1

    # Draw Genes
    gene_count = 0
    for gene in Complete_Gene_List:
        # Draw Gene Line
        gene_height = 100 + (gene_count * 150)
        ct.set_source_rgb(0.1, 0.1, 0.1)
        ct.set_line_width(3)
        ct.move_to(left_border, gene_height)
        ct.line_to(left_border + gene.mLength, gene_height)
        ct.stroke()

        # Draw Exon Box
        ct.set_source_rgb(0.1, 0.1, 0.1)
        ct.set_line_width(3)
        exon_width = gene.mExonStop - gene.mExonStart
        ct.rectangle(left_border + gene.mExonStop, gene_height - 20, exon_width, 40)
        ct.stroke()

        # Draw Gene Header
        header_height = 45 + (gene_count * 150)
        ct.set_source_rgb(0.1, 0.1, 0.1)
        ct.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        ct.set_font_size(15)
        ct.move_to(left_border, header_height)
        ct.show_text(gene.mHeaderName)

        # Draw Motif Lines
        for mt in gene.mMotifList:
            motif_string = temp_dict[mt.mColor]
            motif_color = color_dict[motif_string]
            ct.set_source_rgb(motif_color[0], motif_color[1], motif_color[2])
            ct.set_line_width(3)
            ct.move_to(left_border + mt.mStart_pos, gene_height - 10)
            ct.line_to(left_border + mt.mStart_pos, gene_height  + 10)
            ct.stroke()

        gene_count = gene_count + 1
    surface.write_to_png(out_file)  # Output to PNG

def main():
    longest_gene_length = 0 # Keep track of the longest gene

    multi2linefasta(input_file) # Converts a multi line fasta to a single lined fasta

    f_SaveMotifs(motif_file) # Read motif file and stores all possible motifs into a dictionary

    opened_input = open(input_file, 'r') # open input file

    gene_list = [] # store the genes
    gene_number = 0 # gene counter
    for record in SeqIO.parse(opened_input, "fasta"):
        gene_number = gene_number + 1 # Gene number is one based
        capital_list = []  # list to store index of capital nucleotides
        gene_name = record.id
        header = record.description
        sequence = record.seq
        gene_length = len(sequence)
        if gene_length > longest_gene_length: # Save the longest gene length
            longest_gene_length = gene_length

        for i in range(0, len(sequence)):  # For loop keeps track of the index for every capital letter
            char = sequence[i]
            if char.isupper():
                capital_list.append(i + 1)
        if len(capital_list) == 0:
            exon_start = 0
            exon_stop = 0
        else:
            exon_start = capital_list[0]  # Exon start position is one based
            exon_stop = capital_list[-1]  # Exon end position is one based

        dict_count, list_motif = f_ParseSequence(sequence)

        gene_list.append(cGene(gene_length, exon_start, exon_stop, header, list_motif, gene_number, dict_count))

    f_draw(gene_list, longest_gene_length)

    print("Gene header name followed by the count of motifs in each gene")
    for obj in gene_list:
        print(obj.mHeaderName)
        print(obj.mMotifDict)

    remove_temp_file = "rm singleline.fasta" # Command to move temp file
    os.system(remove_temp_file)  # Remove Temp File

    opened_input.close()  # Close File


main()


