__author__ = 'Kozel'
import re

#fileName = input("Input file path:")
fileName = "/Users/Kozel/Documents/UJ/Papiery Igi/treść.txt"
f = open(fileName, 'r')
string = f.read()
f.close()

# find abbreviations
# assumption: at least 2 capital letters  or 1 capital and 1 number
# no white spaces
# \S nonwhitespace, \s whitespace, [A-Z] - capital, [0-9] numeric
# \w alphanumeric
pattern = r'\s(?P<abb>\w*[A-Z]\S*[A-Z0-9]\w*)\s'

outFile = open('abbreviations2.txt', 'w')

#print(string)
#print(re.findall(pattern, string))
abbreviations = []

for match in re.finditer(pattern, string):
    s = match.group('abb').start()
    e = match.group('abb').end()
    #print(string[s:e])
    if string[s:e] not in abbreviations:
        abbreviations.append(string[s:e])

abbreviations.sort()
#print(abbreviations)
for item in abbreviations:
    outFile.write(item)
    outFile.write('\u2028') #soft return
outFile.close()