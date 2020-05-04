import re

#k => k-mer
def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1): #Initializing
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
    for i in range(n-k+1): #Adding up the frequency
        Pattern = Text[i:i+k]
        freq[Pattern] = freq[Pattern] + 1
    return freq

#Input: A string Text and an integer k
#Output: A list containing all most frequent k-mers in Text
def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if(freq[key] == m):
            words.append(key)
    return words

#OriC of Vibrio Cholero
Text = """ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC"""
k = 9
print(FrequentWords(Text, k))

print(FrequentWords("CGGAGGACTCTAGGTAACGCTTATCAGGTCCATAGGACATTCA", 3))


#--------------Check these to improve------------
#Returns starting position of Pattern in Genome as a list
def PatternMatchingIndex(Pattern, Genome):
    positions = []
    for i in range(len(Genome) - len(Pattern) + 1):
        if Genome[i:i+len(Pattern)] == Pattern:
            positions.append(i)
    return positions

#Pattern counting
def PatternCount(Pattern, String):
    return len(re.findall("(?=%s)"%Pattern, String))


# Regex method
# print the positions variable
# print ([m.start() for m in re.finditer('(?='+ "CTTGATCAT" +')', v_cholerae)])

#Regex to count the no. of occurences
# Finally, print the sum of count_1 and count_2
# print (len(re.findall('(?=' + "ATGATCAAG" + ')', Text)) + len(re.findall('(?=' + "CTTGATCAT" + ')', Text)))