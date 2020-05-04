#Input: A String Genome
#Output: THe skew array of Genome as a list. 
def compute_skew_array(genome):
    skew_array = [0]
    values = {'A':0, 'T':0, 'C':-1, 'G':1}
    for i in range(len(genome)):
        skew_array.append(skew_array[i] + values[genome[i]])
    return skew_array

def compute_minimum_skew(genome):
    skew_array = compute_skew_array(genome)
    min_in_skey_array = min(skew_array)
    positions = []
    for i in range(len(genome)):
        if min_in_skey_array == skew_array[i]:
            positions.append(i)
    return positions

#We say that position i in k-mers p and q is a mismatch if the symbols at position i
# of the two strings are not the same. The total number of mismatches between strings p and q
# is called the Hamming distance between these strings. 
def compute_hamming_distance(p,q):
    num_of_mismatches = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            num_of_mismatches += 1
    return num_of_mismatches

#We say that a k-mer Pattern appears as a substring of Text with at most d mismatches 
# f there is some k-mer substring Pattern' of Text having d or fewer mismatches with Pattern; 
# that is, HammingDistance(Pattern, Pattern') ≤ d. Our observation that a DnaA box may appear 
# with slight variations leads to the following generalization of the Pattern Matching Problem.
def compute_approximate_pattern_matching(Text, Pattern, d):
    positions = [] #Initializing list of positions
    for i in range(len(Text) - len(Pattern) + 1):
        if compute_hamming_distance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions

#This function computes the number of occurrences of Pattern in Text with at most d mismatches.
def compute_approximate_pattern_count(Text, Pattern, d):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if compute_hamming_distance(Text[i:i+len(Pattern)], Pattern) <= d:
            count += 1
    return count


text = 'AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT'
print(compute_skew_array(text))
print(compute_minimum_skey(text))