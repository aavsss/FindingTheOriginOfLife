import random

# Input:
#     AACGTA
#     CCCGTT
#     CACCTT
#     GGATTA
#     TTCCGG
# Output:
#     {'A': [2, 3, 2, 1, 1, 3],
#      'C': [3, 2, 5, 3, 1, 1], 
#      'T': [2, 2, 1, 2, 5, 3], 
#      'G': [2, 2, 1, 3, 2, 2]}
#Adding 1 to each element in the count matrix
#Count-Matrix : No. of nucleotides in each column
#Pseudo-Count-Matrix: No. of nucleodies plus 1 in each column.
#                     Laplace's Rule of Succession: The proabability that sun may never rise tomorrow
def get_pseudo_count_matrix(motifs):
    #Initiaze count matrix as a dictionary
    count = {}

    #Ranging over all nucleotides symbol and creating a list of zeroes
    #corresponding to count[symbol]
    k = len(motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)

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
#     {'A': [0.2222222222222222, 0.3333333333333333, 0.2222222222222222, 0.1111111111111111, 0.1111111111111111, 0.3333333333333333], 
#      'C': [0.3333333333333333, 0.2222222222222222, 0.5555555555555556, 0.3333333333333333, 0.1111111111111111, 0.1111111111111111], 
#      'T': [0.2222222222222222, 0.2222222222222222, 0.1111111111111111, 0.2222222222222222, 0.5555555555555556, 0.3333333333333333], 
#      'G': [0.2222222222222222, 0.2222222222222222, 0.1111111111111111, 0.3333333333333333, 0.2222222222222222, 0.2222222222222222]}
#Dividing each element in count matrix by (no. of items in each column + 4)
#We are adding an extra 1 count to each group in the output (GCAT), not to each element in the motif.  
def get_pseudo_profile_matrix(motifs):
    t = len(motifs)
    k = len(motifs[0])

    psuedo_profile = get_pseudo_count_matrix(motifs)

    for key in psuedo_profile:
        for i in range(k):
            psuedo_profile[key][i] /= (t+4) #t+4 because adding extra 1 to each group in the output

    return psuedo_profile

#t = len(dna)
#k: k-mer
#dna: list of strings dna
#RandomMotifs
#Generates random motifs (k-mers) from each DNA string
# Input:
#     AACGTA
#     CCCGTT
#     CACCTT
#     GGATTA
#     TTCCGG
# Output:
#     AA
#     CG
#     CT
#     GA
#     GG
# Assuming k = 2
def generate_random_motifs(dna,k,t):
    random_kmers = []
    for i in range(t):
        r = random.randint(0,len(dna[0])-k)
        random_kmers.append(dna[i][r:r+k])
    return random_kmers

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
    count_matrix = get_pseudo_count_matrix(motifs)

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

#This function takes a dictionary Probabilities whose keys are k-mers and 
# whose values are the probabilities of these k-mers (which do not necessarily 
# sum to 1). The function should divide each value in Probabilities by the sum
# of all values in  Probabilities, then return the resulting dictionary.
# To make the sum of all probabilities equal to 1
def normalize(probabilities):
    sum_of_probabilites = 0
    for key in probabilities:
        sum_of_probabilites += probabilities[key]
    for key in probabilities:
        probabilities[key] = probabilities[key]/sum_of_probabilites
    return probabilities

#Probabilites are lined up as if they are in a straight line. 
# A random num is generated and is compared against this straight line. 
# The corresponding key is returned. 
# In-short: A random motif is selected based on a weighted probability scale.  
def weighted_die(probabilities):
    random_num = random.uniform(0,1)
    count = 0
    for key in probabilities:
        count += probabilities[key]
        if random_num < count:
            return key

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

#Assemble this code into a function ProfileGeneratedString(Text, profile, k) that takes a string Text,
#  a profile matrix profile , and an integer k as input.  
# It should then return a randomly generated k-mer from Text whose probabilities are generated from profile, as described above.  
# Computes a motif (k-mer) with the most probability based on a given proabability and a DNA String. 
#Sample Input:
# 1
# AAACCCAAACCC
# {'A': [0.5, 0.1], 
#  'C': [0.3, 0.2], 
#  'G': [0.2, 0.4], 
#  'T': [0.0, 0.3]}
# 2
# Sample Output:
# AC
def profile_generated_string(Text, profile_matrix,k):
    n = len(Text)
    probabilities = {}
    for i in range(n-k+1):
        probabilities[Text[i:i+k]] = compute_probability(Text[i:i+k], profile_matrix)
    probabilities = normalize(probabilities)
    return weighted_die(probabilities)

# Generates a random motif at first
# It is assigned as the best motif
# In each iteration:
# Delete a random row of motif.
# Calculate pseudo profile matrix from the remaining motif
# Generate the most probable motif based on weighted probability and replace the deleted motif
# Compare the score, if the new motifs scores less then it is assigned as the best motif
# If there were 10 rows then there would be 10 iterations
def gibbs_sampler(dna,k,t,n):
    best_motifs = []
    motifs = generate_random_motifs(dna,k,t)
    best_motifs = motifs
    for j in range(n):
        i = random.randint(0,t-1)
        del motifs[i]
        profile_matrix = get_pseudo_profile_matrix(motifs)
        motifs.insert(i, profile_generated_string(dna[i], profile_matrix, k))
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

#setting up variables:
#t : Number of strings in dna
#k : k-mer
#N : Number of times we want to do gibbs sampling
t = 10
k = 15
N = 100

#Calling gibbs sampling 100 times
#Storing the least scoring (best) motif
def repeated_gibbs_sampler(dna,k,t,N):
    best_score = float('inf')
    best_motifs = []
    for i in range(N):
        motifs = gibbs_sampler(dna,k,t,N)
        curr_score = compute_score(motifs)
        if curr_score < best_score:
            best_score = curr_score
            best_motifs = motifs
    return best_motifs

best_motifs = repeated_gibbs_sampler(dna,k,t,N)
print(best_motifs)
print(compute_score(best_motifs))
