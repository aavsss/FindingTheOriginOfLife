import random

#Pre-condition: A list of strings as sent as an input
# Input:
#     AACGTA
#     CCCGTT
#     CACCTT
#     GGATTA
#     TTCCGG
# Output:
#     {'A': [1, 2, 1, 0, 0, 2], 
#      'C': [2, 1, 4, 2, 0, 0], 
#      'G': [1, 1, 0, 2, 1, 1], 
#      'T': [1, 1, 0, 1, 4, 2]}
# Computes number of nucleotides in each column
def get_count_matrix(motifs):
    #Initiaze count matrix as a dictionary
    count = {}

    #Ranging over all nucleotides symbol and creating a list of zeroes
    #corresponding to count[symbol]
    k = len(motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)

    #Ranging over all elements symbol=Motifs[i][j] of the count matrix 
    #and add 1 to count[symbol][j]
    t = len(motifs)
    for i in range(t):
        for j in range(k):
            symbol = motifs[i][j]
            count[symbol][j] += 1
    return count

# Input:
#     AACGTA
#     CCCGTT
#     CACCTT
#     GGATTA
#     TTCCGG
# Output:
#     {'A': [0.2, 0.4, 0.2, 0.0, 0.0, 0.4], 
#      'C': [0.4, 0.2, 0.8, 0.4, 0.0, 0.0], 
#      'G': [0.2, 0.2, 0.0, 0.4, 0.2, 0.2], 
#      'T': [0.2, 0.2, 0.0, 0.2, 0.8, 0.4]}
# Computed by divided the number of nucleotides in each column by the total no. of rows. 
# Sum of each column should be 1
def get_profile_matrix(motifs):
    t = len(motifs)
    k = len(motifs[0])

    profile = get_count_matrix(motifs)

    for key in profile:
        for i in range(k):
            profile[key][i] = profile[key][i]/t

    return profile


# Takes the count matrix and prints out a string which has the most value in each column
# If two nucleotides have same count number then a nucleotide is selected randomly. 
# Input:
#     AACGTA
#     CCCGTT
#     CACCTT
#     GGATTA
#     TTCCGG
# Count Matrix:
#     {'A': [2, 3, 2, 1, 1, 3],
#      'C': [3, 2, 5, 3, 1, 1], 
#      'T': [2, 2, 1, 2, 5, 3], 
#      'G': [2, 2, 1, 3, 2, 2]}
# Output:
#  Consensus String:
#   CACCTA
def get_consensus_string(motifs):
    k = len(motifs[0])
    count_matrix = get_count_matrix(motifs)

    #The string
    consensus = ""
    for j in range(k):
        m = 0
        frequent_symbol = ""
        for symbol in "ACGT":
            if count_matrix[symbol][j] > m:
                m = count_matrix[symbol][j]
                frequent_symbol = symbol
        consensus += frequent_symbol

    return consensus

# Takes in the consensus string and computes how many nucleotides are off in each column
# Input:
#     AACGTA
#     CCCGTT
#     CACCTT
#     GGATTA
#     TTCCGG
# Consensus String: CACCTA 
# Output:
#     3+3+1+3+1+3  
#     14
def compute_score(motifs):
    consensus_string = get_consensus_string(motifs)
    k = len(motifs)
    l = len(motifs[0])
    score = 0

    for i in range(k):
        for j in range(l):
            if motifs[i][j] != consensus_string[j]:
                score += 1
    
    return score

# Input:
#     ACGGGGATTACC
#     0.2 0.2 0.0 0.0 0.0 0.0 0.9 0.1 0.1 0.1 0.3 0.0
#     0.1 0.6 0.0 0.0 0.0 0.0 0.0 0.4 0.1 0.2 0.4 0.6
#     0.0 0.0 1.0 1.0 0.9 0.9 0.1 0.0 0.0 0.0 0.0 0.0
#     0.7 0.2 0.0 0.0 0.1 0.1 0.0 0.5 0.8 0.7 0.3 0.4
# Output:
#     0.000839808
# To implement a function Pr(Text, Profile), we begin by 
# setting a “probability” variable p equal to 1. We then range
# through the characters of Text one at a time. At position i of Text, 
# we set p equal to p times the value of Profile corresponding to symbol
# Text[i] and column i, which is just Profile[Text[i]][i].
#The probability that a profile matrix will produce a given string is given by the product of individual nucleotide probabilities.
#Computes probability of the k-mer occuring. 
# Multiplies the probability of nucleotide in each column
def compute_probability(text, profile_matrix):
    probability = 1
    for i in range(len(text)):
        probability *= profile_matrix[text[i]][i]
    return probability

#Computes the most probable kmer based on the profile matrix
# Compares probability of each kmer
def compute_profile_most_probable_kmer(text, k, profile_matrix):
    max_probability = -1
    result = ""
    for i in range(len(text)-k+1):
        k_mer = text[i:i+k]
        temp_probability = compute_probability(k_mer, profile_matrix)
        if temp_probability > max_probability:
            result = k_mer
            max_probability = temp_probability
    return result

# Test 0 # Sample Dataset 
# Input:
#     3 5
#     GGCGTTCAGGCA
#     AAGAATCAGTCA
#     CAAGGAGTTCGC
#     CACGTCAATCAC
#     CAATAATATTCG
# Output:
#     CAG
#     CAG
#     CAA
#     CAA
#     CAA
# t = number of k-mers in dna
#http://www.mrgraeme.co.uk/greedy-motif-search/
#Steps:
# First k-mer is extracted from each string and is assigned as the best motif
# In each iteration of seletecting every kmer from the string
# New profile matrix is build after making motifs from each string(row) of DNA
# Next motif is added based on the newly created profile matrix
# Motif with the least score is iniialized as the best motif. 
def search_greedy_motif(dna,k,t):
    best_motifs = []
    for i in range(1,t):
        best_motifs.append(dna[i][0:k])
    
    n = len(dna[0])
    for i in range(n-k+1):
        motifs = []
        motifs.append(dna[0][i:i+k])
        for j in range(1,t):
            p = get_profile_matrix(motifs[0:j])
            motifs.append(compute_profile_most_probable_kmer(dna[j], k, p))
        if compute_score(motifs) < compute_score(best_motifs):
            best_motifs = motifs
    return best_motifs

#Hyperlinked DosR dataset
dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC",
        "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG",
        "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC",
        "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC",
        "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG",
        "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA",
        "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA",
        "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG",
        "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG",
        "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]

t = len(dna)
k = 15
motifs = search_greedy_motif(dna,k,t)
print(motifs)
print(compute_score(motifs))