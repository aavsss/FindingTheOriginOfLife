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
#We are adding an extra 1 count to each group in the output (GCAT), not to each element in the motif.  
def get_pseudo_profile_matrix(motifs):
    t = len(motifs)
    k = len(motifs[0])

    psuedo_profile = get_pseudo_count_matrix(motifs)

    for key in psuedo_profile:
        for i in range(k):
            psuedo_profile[key][i] /= (t+4) #t+4 because adding extra 1 to each group in the output

    return psuedo_profile

# Input:
#     AACGTA
#     CCCGTT
#     CACCTT
#     GGATTA
#     TTCCGG
# Output:
#     CACCTA

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

    # Input:
#     AACGTA
#     CCCGTT
#     CACCTT
#     GGATTA
#     TTCCGG
# Output:
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
def compute_probability(text, profile_matrix):
    probability = 1
    for i in range(len(text)):
        probability *= profile_matrix[text[i]][i]
    return probability

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

# Test 0 # Sample Dataset (your code is not run on this dataset)
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
def search_greedy_motif(dna,k,t):
    best_motifs = []
    for i in range(1,t):
        best_motifs.append(dna[i][0:k])
    
    n = len(dna[0])
    for i in range(n-k+1):
        motifs = []
        motifs.append(dna[0][i:i+k])
        for j in range(1,t):
            p = get_pseudo_profile_matrix(motifs[0:j])
            motifs.append(compute_profile_most_probable_kmer(dna[j], k, p))
        if compute_score(motifs) < compute_score(best_motifs):
            best_motifs = motifs
    return best_motifs

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
def compute_motifs(profile_matrix, dna, k):
    motifs = []
    for i in range(len(dna)):
        motifs.append(compute_profile_most_probable_kmer(dna[i], k, profile_matrix))
    return motifs

#t = len(dna)
#k: k-mer
#dna: list of strings dna
#RandomMotifs
def generate_random_kmer_from_strings(dna,k,t):
    random_kmers = []
    for i in range(t):
        r = random.randint(0,len(dna[0])-k)
        random_kmers.append(dna[i][r:r+k])
    return random_kmers

def compute_randomized_motif_search(dna,k,t):
    m = generate_random_kmer_from_strings(dna,k,t)
    best_motifs = m

    while True:
        profile = get_pseudo_profile_matrix(m)
        m = compute_motifs(profile,dna)
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

# BestMotifs = RepeatedRandomizedMotifSearch(dna, k, t)
# # Print the BestMotifs variable
# print(BestMotifs)
# print(compute_score(BestMotifs))

