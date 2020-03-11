# Motif Mark

***Author:*** Wyatt Eng

***Goal:*** Parse a FASTA file that contains a gene (with introns and exons) and identify biological motifs. Motifs are short sequences of DNA that have biological significance. The current program is scalable and can handle any number of sequences and motifs. 

### Required Input
- -f : The path to a FASTA file containing genes/mRNAs, with exons in uppercase and introns in lowercase.
- -m : A text file with a list of motifs. Each motif should be on an individual line. Motifs may contain IUPAC ambiguous nucleotides.
- -o : The path and name of the output png file.
- -t : Molecule type of input fasta file. DNA (contains Thymine) or RNA (contains Uracil).
