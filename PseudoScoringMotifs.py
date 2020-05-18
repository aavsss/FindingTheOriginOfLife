import random

# Input:
#     AACGTA
#     CCCGTT
#     CACCTT
#     GGATTA
#     TTCCGG
# Output:
#     {'A': [2, 3, 2, 1, 1, 3], 'C': [3, 2, 5, 3, 1, 1], 'T': [2, 2, 1, 2, 5, 3], 'G': [2, 2, 1, 3, 2, 2]}
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
#     {'A': [0.2222222222222222, 0.3333333333333333, 0.2222222222222222, 0.1111111111111111, 0.1111111111111111, 0.3333333333333333], 'C': [0.3333333333333333, 0.2222222222222222, 0.5555555555555556, 0.3333333333333333, 0.1111111111111111, 0.1111111111111111], 'T': [0.2222222222222222, 0.2222222222222222, 0.1111111111111111, 0.2222222222222222, 0.5555555555555556, 0.3333333333333333], 'G': [0.2222222222222222, 0.2222222222222222, 0.1111111111111111, 0.3333333333333333, 0.2222222222222222, 0.2222222222222222]}
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

# NEXT SECTION
# Input:
#     0.8 0.0 0.0 0.2
#     0.0 0.6 0.2 0.0
#     0.2 0.2 0.8 0.0
#     0.0 0.2 0.0 0.8
#     TTACCTTAAC
#     GATGTCTGTC
#     ACGGCGTTAG
#     CCCTAACGAG
#     CGTCAGAGGT
# Output:
#     ACCT
#     ATGT
#     GCGT
#     ACGA
#     AGGT
# Returns most probable motifs based on the profile matrix
def compute_motifs(profile_matrix, dna, k):
    motifs = []
    for i in range(len(dna)):
        motifs.append(compute_profile_most_probable_kmer(dna[i], k, profile_matrix))
    return motifs

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

#Randomly generate motifs from strings of DNA
# Assign it as the best motifs
# In iteration:
#  Iteration is done until the value of score is minimized. 
#  Generate a profile matrix based on the motifs computed
#  motifs is computed from the most probable kmer based on the profile matrix in each DNA String
#  If the score improves (becomes less) then the iteration continues
#  Else the best motif computed is returned. 
def compute_randomized_motif_search(dna,k,t):
    m = generate_random_motifs(dna,k,t)
    best_motifs = m

    while True:
        profile = get_pseudo_profile_matrix(m)
        m = compute_motifs(profile,dna,k)
        if compute_score(m) < compute_score(best_motifs):
            best_motifs = m
        else:
            return best_motifs

#Call RandomizedMotifSearch(dna,k,t) N times, storing the best-scoring set of motifs
def RepeatedRandomizedMotifSearch(dna, k, t):
    BestScore = float('inf')
    BestMotifs = []
    N = 100
    for i in range(N):
        Motifs = compute_randomized_motif_search(dna, k, t)
        CurrScore = compute_score(Motifs)
        if CurrScore < BestScore:
            BestScore = CurrScore
            BestMotifs = Motifs
    return BestMotifs

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

BestMotifs = RepeatedRandomizedMotifSearch(dna, k, t)
# Print the BestMotifs variable
print(BestMotifs)
print(compute_score(BestMotifs))


