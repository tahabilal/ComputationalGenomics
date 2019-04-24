
from random import randint
import numpy as np 

input_text = input("Enter a file: ")
if len(input_text) < 1 : input_text = "dna_sample.txt"
Text = open(input_text)
Text = Text.read().split()


k = int(input("Enter k: "))

def random_motif(a):
    motifs = []
    for i in range(0,10):
        y = randint(0, 500-a)
        motifs.append(Text[i][y:y+a])
    return motifs


def randomized_motif_search(randomMotif,a): 
    A, C, G, T = [], [], [], []
    consensus = ""
    transpose_motif = [" "," "," "," "," "," "," "," "," "," "," "]
    for j in range(a):
        countA, countC, countG, countT = 0, 0, 0, 0
        for i in range(10):
            transpose_motif[j] += randomMotif[i][j]
            if randomMotif[i][j] == "A" or randomMotif[i][j] == "a":
                countA = countA + 1
            if randomMotif[i][j] == "C" or randomMotif[i][j] == "c":
                countC = countC + 1
            if randomMotif[i][j] == "G" or randomMotif[i][j] == "g":
                countG = countG + 1
            if randomMotif[i][j] == "T" or randomMotif[i][j] == "t":
                countT = countT + 1
        if countA >= max(countC, countG, countT):
            consensus = consensus + "A"
        elif countC >= max(countA, countG, countT):
            consensus = consensus + "C"
        elif countG >= max(countA, countC, countT):
            consensus = consensus + "G"
        elif countT >= max(countA, countC, countG):
            consensus = consensus + "T"
    
        A.append(countA)
        C.append(countC)
        G.append(countG)
        T.append(countT)

    count_matrix = np.array([ A,C,G,T ])
    score = calculate_score(transpose_motif, consensus,a)
    profile_matrix = np.true_divide(count_matrix, len(A))
    updated_motifs = calculate_new_motifs(profile_matrix, Text, a)
    return count_matrix, score, profile_matrix, updated_motifs, consensus


def calculate_score(transpose_motif, consensus,a):
    score = 0
    for j in range(10):
        for i in range(a):
            if transpose_motif[j][i].upper() != consensus[j]:
                score += 1
    return score


def calculate_new_motifs(profile_matrix, Text,a):
    patterns = []
    patterns_score = []
    pattern_score = 1
    updated_motifs = []
    index = 0
    for j in range(10):
    
        for i in range(len(Text[j])-a):
            if i+a > (500):
                break
            pattern = Text[j][i: i+a]
            patterns.append(pattern)
    
        for i in range(len(patterns)):
            for x in range(len(patterns[i])):
                if patterns[i][x].upper() == "A":
                    pattern_score *= profile_matrix[0][x]
                if patterns[i][x].upper() == "C":
                    pattern_score *= profile_matrix[1][x]
                if patterns[i][x].upper() == "G":
                    pattern_score *= profile_matrix[2][x]
                if patterns[i][x].upper() == "T":
                    pattern_score *= profile_matrix[3][x]
            patterns_score.append(pattern_score)
            pattern_score = 1
    
        max_score = max(patterns_score)
    
        for y in range(len(patterns_score)):
            if patterns_score[y] == max_score:
                index = y
                break
    
        updated_motifs.append(patterns[index])
        patterns_score = []
        patterns = []
    
    return updated_motifs


def calculate_new_motifs_gibbs(profile_matrix, Text, index, randomMotif,a):  # new method for gibbs sampler
    patterns = []
    patterns_score = []
    pattern_score = 1
    updated_motifs = randomMotif
    for i in range(len(Text[index])):
        if i + a > 500:
            break
        pattern = Text[index][i: i + a]
        patterns.append(pattern)

    for i in range(len(patterns)):
        for y in range(len(patterns[i])):
            if patterns[i][y].upper() == "A":
                pattern_score *= profile_matrix[0][y]
            if patterns[i][y].upper() == "C":
                pattern_score *= profile_matrix[1][y]
            if patterns[i][y].upper() == "G":
                pattern_score *= profile_matrix[2][y]
            if patterns[i][y].upper() == "T":
                pattern_score *= profile_matrix[3][y]
        patterns_score.append(pattern_score)
        pattern_score = 1

    # probabilities' sum must be 1
    patterns_score = np.array(patterns_score)
    patterns_score /= patterns_score.sum()
    updated_motifs[index] = np.random.choice(patterns, p=patterns_score)  # random choose according to profile matrix

    return updated_motifs


def construct_consensus(randomMotif,a):  # new method for gibbs sampler
    consensus = ""
    transpose_motif = [" "," "," "," "," "," "," "," "," "," "," "]
    for j in range(a):
        countA, countC, countG, countT = 0, 0, 0, 0
        for i in range(10):
            transpose_motif[j] +=randomMotif[i][j] 
            if randomMotif[i][j] == "A" or randomMotif[i][j] == "a":
                countA = countA + 1
            if randomMotif[i][j] == "C" or randomMotif[i][j] == "c":
                countC = countC + 1
            if randomMotif[i][j] == "G" or randomMotif[i][j] == "g":
                countG = countG + 1
            if randomMotif[i][j] == "T" or randomMotif[i][j] == "t":
                countT = countT + 1
        if countA >= max(countC, countG, countT):
            consensus = consensus + "A"
        elif countC >= max(countA, countG, countT):
            consensus = consensus + "C"
        elif countG >= max(countA, countC, countT):
            consensus = consensus + "G"
        elif countT >= max(countA, countC, countG):
            consensus = consensus + "T"
    return consensus, transpose_motif


def gibbs_sampler(randomMotif,a):  # new method for gibbs sampler
    A, C, G, T = [], [], [], []
    index = randint(0, 9)
    consensus, transpose_motif = construct_consensus(randomMotif,a)
    print(consensus)
    score = calculate_score(transpose_motif, consensus,a)

    for j in range(a):
        countA, countC, countG, countT = 0, 0, 0, 0
        for i in range(10):
            if i == index:  # remove one motif as a random (check line 165)
                continue
            if randomMotif[i][j] == "A" or randomMotif[i][j] == "a":
                countA = countA + 1
            if randomMotif[i][j] == "C" or randomMotif[i][j] == "c":
                countC = countC + 1
            if randomMotif[i][j] == "G" or randomMotif[i][j] == "g":
                countG = countG + 1
            if randomMotif[i][j] == "T" or randomMotif[i][j] == "t":
                countT = countT + 1
        A.append(countA+1)
        C.append(countC+1)
        G.append(countG+1)
        T.append(countT+1)

    count_matrix = np.array([A, C, G, T])
    profile_matrix = np.true_divide(count_matrix, len(A))
    updated_motifs = calculate_new_motifs_gibbs(profile_matrix, Text, index, randomMotif,a)
    return count_matrix, score, profile_matrix, updated_motifs, consensus


def exp1():
    best_score = 10000
    iteration = 0
    randomMotifs = random_motif(k)
    print(randomMotifs)
    
    while 1:
        count_matrix, score, profile_matrix, updated_motifs, consensus = randomized_motif_search(randomMotifs,k)
        randomMotifs = updated_motifs
        if score < best_score:
            best_score = score
            print(profile_matrix)
            print("Consensus String:",consensus)
            print("Score:", score)
            final_motifs = updated_motifs
            final_consensus = consensus
            final_score = score
            iteration += 1
        else:
            print("\n---- Results ----")
            print("Algorithm executed in {} iterations".format(iteration))
            print("\nFinal Motifs:")
            print(*final_motifs, sep='\n')
            print("\nConsensus String:")
            print(final_consensus)
            print("\nFinal Score:")
            print(final_score)
            break


def exp2():  # new method for gibbs sampler
    randomMotifs = random_motif(k)
    print(randomMotifs)
    count_non_improve = 0
    best_score = 10000

    while True:
        count_matrix, score, profile_matrix, updated_motifs, consensus = gibbs_sampler(randomMotifs,k)
        randomMotifs = updated_motifs
        if score < best_score:
            best_score = score
            print(profile_matrix)
            print("Consensus String:", consensus)
            print("Score:", score)
            print("-----------------")
            print("Updated Motifs:")
            print("-----------------")
            print(*updated_motifs, sep='\n')
            print("\n")
            final_motifs = updated_motifs
            final_consensus = consensus
            final_score = score
            count_non_improve = 0
        else:
            print("New Score calculated as: {} ".format(score))
            print("The Score didn't improved since last {}".format(count_non_improve))
            count_non_improve = count_non_improve +1
            if count_non_improve == 50:
                break
    print("\n---- Results ----")
    print("\nFinal Motifs:")
    print(*final_motifs, sep='\n')
    print("\nConsensus String:")
    print(final_consensus)
    print("\nFinal Score:")
    print(final_score)
            
            

def main():
    print("Menu\n1-Randomized Motif Search\n2-Gibbs Sampler\n")
    case = input("Enter your choice: ")
    while 1:
        if case == "1":
            exp1()
        elif case == "2":
            exp2()
        else:
            print("Terminated")
            exit(0)
        case = input("\nEnter your choice: ")


if __name__ == "__main__":
    main()


   
























