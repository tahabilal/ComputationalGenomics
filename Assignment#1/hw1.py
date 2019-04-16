

from collections import Counter
import re

print("\nDefault file name(InputA) is assigned when it was not entered any file name\n")
input_text = input("Enter a file name (without \".txt\"): ") + ".txt"
if len(input_text) < 1:
    input_text = "InputA.txt"
Text = open(input_text)
Text = Text.read()
Text = Text.upper()

k = int(input("Input k: "))
x = int(input("Input x: "))

# error check  0 < k < 10    2 < x
while k < 1 or k > 9 or x < 2:
    print("Please enter greater than 1 for k and greater than 2 for x or less than 9 for only k")
    k = int(input("Input k: "))
    x = int(input("Input x: "))


def frequent_words_problem(Text,k):

        patterns = []  # all patterns in there
        filter_patterns = []  # all greater than x value patterns
        for i in range(len(Text)):
            pattern = Text[i: i + k]
            if len(pattern) == k:
                patterns.append(pattern)
        count = Counter(patterns).most_common()
        for i in range(len(count)):
            if count[i][1] >= x:  # counts are greater than or equal to x value.
                filter_patterns.append(count[i][0])
        filter_patterns_str = ",".join(filter_patterns)
        return filter_patterns, filter_patterns_str


def reverse_complement(pattern):
    # find reverse complement given pattern
    pattern = pattern[::-1] # reverse string
    revComp = {"A":"T", "T":"A", "G":"C", "C":"G"}
    newlist = []
    for c in pattern:
        newlist.append(revComp[c])
    return "".join(newlist)


freq_words = frequent_words_problem(Text, k)
complement = [reverse_complement(v) for v in freq_words[0]]
occurrences = []
for c in complement:
    tmp = [m.start() for m in re.finditer('(?={})'.format(c),Text)]
    occurrences.append(len(tmp))
print("Outputs:\n"+ str(k) + "-mer: " + freq_words[1])
print("Reverse complement: ")
for i, j in zip(complement,occurrences):
    print(i + " appearing " + str(j) + " times")




