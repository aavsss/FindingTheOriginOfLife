
#Input: A string Pattern
#Output: The reverse of Pattern
def Reverse(Pattern):
    rev = ""
    for char in Pattern:
        rev = char + rev
    return rev

def ReversePython(Pattern):
    return Pattern[::-1]

# Input:  A DNA string Pattern
compdict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
# Output: The complementary string of Pattern (with every nucleotide replaced by its complement).
def Complement(Pattern):
    # your code here
    compstring = ''
    for char in Pattern:
        compstring = compstring + compdict[char]
    return compstring

def ComplementAlternate(Pattern):
    return ''.join([compdict[n] for n in Pattern])

#Reverse Complement function because we read 5' to 3'
def ReverseComplement(Pattern):
    Pattern = Reverse(Pattern)
    Pattern = Complement(Pattern)
    return Pattern

print(ReverseComplement("GATTACA"))