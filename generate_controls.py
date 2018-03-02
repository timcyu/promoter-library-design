# This script generates positive control promoters that will be added to the library
import re
import argparse
import random
import subprocess


# This function extracts the promoters from a CSV file
def extract_syn_sequences(filename):
    infile = open(filename, 'r')
    syn_promoters = {}

    # read through header
    infile.readline()

    for line in infile.readlines():
        fields = line.strip().split(',')
        name = fields[0]
        seq = fields[9]

        # remove white space from seq
        seq = ''.join(seq.split())
        # remove quotations
        match = re.search('[ACGT]{1,}', seq)
        clean_seq = match.group(0)

        # remove quotations from name, always first and last characters
        clean_name = name.replace('\"', '')

        syn_promoters['pos_control_' + clean_name] = clean_seq

    return syn_promoters


# This function checks the oligo and its reverse complement and chooses the one with the least amount of A's
# A's are harder to synthesize
def best_A_content(oligo):
    countA = 0
    countA_rev = 0
    for i in oligo:
        if i == 'A':
            countA = countA + 1
    for i in reverse_complement(oligo):
        if i == 'A':
            countA_rev = countA_rev + 1
    if countA_rev > countA:
        final_oligo = oligo
    else:
        final_oligo = reverse_complement(oligo)
    return final_oligo


# All sequences must be the same length. This function adds stuffer sequence to ensure the promoters are the same length
# Forward Primer + Stuffer + XhoI site + Promoter + Reverse Primer = 225 bp
def add_stuffer(syn_promoters, fwd_primer, rev_primer):
    xhoI = 'CTCGAG'
    stuffed_promoters = {}
    stuffer = "CCCGCGCGTTGGCCGATTCATTAATGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGCCCGCGCGTTGGCCGATTCATTAACCCGCGCGTTGGCCGATTCATTAATGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCGTGCAGCTGGCACGACAGGTTTCCCGACTGGAAAGCG"

    for x in syn_promoters:
        promoter = syn_promoters[x]
        to_add = 177 - len(promoter)
        new_promoter = best_A_content(fwd_primer[:20] + stuffer[:to_add - 2] + xhoI + promoter + rev_primer)
        stuffed_promoters[x] = new_promoter
    return stuffed_promoters


# Function returns reverse complement of a string of nucleotides
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 't': 'a', 'a': 't', 'c': 'g', 'g': 'c'}
    rc = ''.join([complement[nt] for nt in seq[::-1]])
    return rc


# Function creates a file containing the final sequences. File format based on argument input
def seq_writer(seqs, output_name, output_type):
    outfile = open(output_name, 'w')

    for x in seqs:

        if output_type == 'fasta':
            outfile.write('>' + x + '\n' + seqs[x] + '\n')
        if output_type == 'tab':
            outfile.write('>' + x + '\t' + seqs[x] + '\n')

    outfile.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate control sequences for library')
    parser.add_argument('syn_file', help='CSV of PNAS synthetic promoters')
    parser.add_argument('fwd_primer', help='FASTA file of forward primer sequence')
    parser.add_argument('rev_primer', help='FASTA file of reverse primer sequence')

    args = parser.parse_args()

    syn_promoters = extract_syn_sequences(args.syn_file)

    with open(args.rev_primer) as infile:
        # read through header
        header = infile.readline()
        rev_primer = infile.readline().strip()
    with open(args.fwd_primer) as infile:
        header = infile.readline()
        fwd_primer = infile.readline().strip()

    stuffed_promoters = add_stuffer(syn_promoters, fwd_primer, rev_primer)

    print "Number of control sequences:", len(stuffed_promoters)
    for x in stuffed_promoters:
        print len(stuffed_promoters[x])

    seq_writer(stuffed_promoters, 'control_sequences.txt', 'tab')
