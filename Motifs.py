#Motifs = ["AACGTA", "CCCGTT", "CACCTT", "GGATTA", "TTCCGG"]


# Input:  A set of kmers Motifs
# Output: Count(Motifs), {'A': [1, 2, 1, 0, 0, 2], 'G': [1, 1, 0, 2, 1, 1], 'T': [1, 1, 0, 1, 4, 2], 'C': [2, 1, 4, 2, 0, 0]}
def Count(Motifs):
    count = {} # initializing the count dictionary
    k = len(Motifs[0])  # sets k equal to the length of Motifs[0], the first string in Motifs
    t = len(Motifs)
    for symbol in "ACGT":
        count[symbol] = []  # empty list
        for j in range(k):
             count[symbol].append(0)

    for i in range(t):      # repeat for t many motifs/strings
        for j in range(k):  # repeat for the length of one motif/string
            symbol = Motifs[i][j]
            count[symbol][j] += 1

    return count
#print("Count returns:", Count(Motifs))

# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus(Motifs):
    k = len(Motifs[0])      # num characters in each string
    count = Count(Motifs)
    consensus = ""
    for j in range(k):      # go through Values list 'A': [1, 2, 1, 0, 0, 2]
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:    # go through 4 letters at Value list position [j]
                m = count[symbol][j]
                frequentSymbol = symbol # assign freq letter
        consensus += frequentSymbol     # update string
    return consensus
#Consensus(Motifs)
#print(Consensus(Motifs))

# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
def Score(Motifs):
    k = len(Motifs[0])      # num characters in each string
    t = len(Motifs)

    consensus = Consensus(Motifs)
    mismatch = 0
    for j in range(k):                      # traverse string
        for i in range(t):                  # traverse list of strings
            symbol = Motifs[i][j]
            if symbol != consensus[j]:
                mismatch += 1
    return mismatch
#print(Score(Motifs))

# Input: dictionary of lists with counts of Motifs
# Output: dictionary of lists with profile matrix of Motifs
def Profile(Motifs):
    t = len(Motifs)         # num strings
    k = len(Motifs[0])      # num characters in each string
    profile = {}            # empty dict
    profile = Count(Motifs)
    # to get a list of value list, use dic.values()
    # dic.items() returns pair list of key and values
    stringList = profile.values()
    for string in stringList:
        for i in range(len(string)):
            string[i] = string[i] / t
    return profile
#print(Profile(Motifs))

# Returns probability of a given sequence (Text) by reading through a probability matrix (Profile)
# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile), probability
def Pr(Text, Profile):
    p = 1.0
    for j in range(len(Text)):
        base = Text[j]              # select one base
        p_base = Profile[base][j]   # find the key corresponding to the base in dict,
                                    #  and a value corresponding to string length
        p *= p_base                 # update p
    return p

#Profile = {'A': [0.273, 0.318, 0.318, 0.333, 0.303, 0.273, 0.227, 0.258, 0.273, 0.167, 0.121, 0.288, 0.227, 0.242, 0.227], 'C': [0.242, 0.212, 0.227, 0.152, 0.197, 0.167, 0.258, 0.227, 0.318, 0.242, 0.273, 0.303, 0.242, 0.197, 0.242], 'G': [0.242, 0.258, 0.242, 0.197, 0.227, 0.303, 0.273, 0.273 ,0.167, 0.258, 0.288 ,0.227 ,0.258 ,0.288, 0.273], 'T': [0.242, 0.212, 0.212, 0.318, 0.273 ,0.258 ,0.242 ,0.242 ,0.242, 0.333, 0.318 ,0.182 ,0.273, 0.273, 0.258]}
#Text = "AGCAGTGGTAGTGATTTGTTTAGTTTGGCATAATTAGAATACCTAGATGTTCGGAGGCCCTCATACTATAGAAAAGGCTGAGATGTCCCTAAGTATGTCAGGTTTGAAAGTAAGGCCATGGGCGAGTCAAGGGTTTACAAGCGTGCATCGGCTCATCGGCTCACAGGTCGTAACTGAGGCTCTGCGAACGGTTACACTCCCCTGGACAATTGCGGTAGTTTTTTGATGAACAGCAAACTGACTGGTACTCTGCCTCCAACGGAATTCTAATCGTGCCTCGACTAGAAGTTAGATCTATTACTTGCCTGGCTGGCACGTTAAAAAGGACTTACCAAGCTAACACGGCATGCCACTACGTGTAAAGCAGAGCACTACCGTAAGCGTTTGTGACGTCTATAAACCGACCATGCTGTTGAGGTGAGCTATCGTAAAACCGCCTATCTGTGGGTTGAAGAAAACACCACTATACGACGTGCCCGCTTTGCTTCGAAGGCTGCTGAGCGGATTAGCTACATAGCTTCACCGCAGGCCCTGATTTGATTTTACGTGTATCTCGTAAGTAAGTTCAGCACGGAGTGTGCTTGACACATCCTTCTCGCAGATGCACAAATGTTCTTTTACCAGGTGTTGAATGCCGTATCACTCGTTAGTGGGGCCAGCACAACTCTTCAATGAAGATGGGAGCTGGAACCAGTACTTATACCCGTGTTTGAGGACGGGCTTTCGGGCGATTCATGAATAACATAGCTGCTTAGGCGAGGCTGGTAAACCTACCATGTCAGAGAAGCTATTACTGCTCACCCCGCTGTGGTACAGCTAAGACAGCCCTAGGCGCATAAAAAGACACTGGTGTCTGTCCTTCGTTAAACCTAGTTTTTTTTGAGGTACAAATTCGTGATGTAGGCTAATATCAAACCTCAGTATCTGGTTTACATATGACAAAACTGACGTGTCGCGGCCTTAGCTGGAATCGAGCTTGCGTTTCTGGGACATCAAAGTA"
#print(Pr(Text,Profile))

# Slides a window of length k across sequence, and evaluates the probability of window seq matching the Profile matrix
# Input:  String Text, an integer k, and profile matrix Profile
# Output: ProfileMostProbablePattern(Text, k, Profile)
def ProfileMostProbablePattern(Text, k, Profile):
    p_max = -1
    for i in range(0, len(Text) - k + 1):
        k_text = Text[i:i+k]            # grab a text of length k
        p_k_text = Pr(k_text, Profile)  # calculate Pr for that short piece
        if p_k_text > p_max:
            p_max = p_k_text
            best_text = k_text
    return best_text
k = 15
#print("here:", ProfileMostProbablePattern(Text, k, Profile))

# https://hahamo.wordpress.com/2015/06/23/greedy-motif-search/
# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    # length of the first Dna string
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))

        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs

    return BestMotifs

#Dna = ["","",""]
Dna = [line.strip() for line in open("test4.txt", 'r')]
print(Dna)

k = 12
t = 25
answer = GreedyMotifSearch(Dna, k, t)
for i in answer: print(i)

# Input: A collection of strings Dna, and integers k and d.
# Output: All (k, d)-motifs in Dna.
def MotifEnumeration(Dna, k, d):
    Patterns = []
