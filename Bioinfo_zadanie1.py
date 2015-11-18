__author__ = 'Kozel'
# coding = utf-8
from Bio.Data import CodonTable
from Bio import SeqIO
import re
'''
część 1
Używając tabeli z kodem genetycznym, przetłumacz następujące wyrażenie
z sekwencji aminokwasowej na sekwencję nukleotydową:
G{1,3}.[^P]{2}[RGI]*L+V
'''
print("\nCzęść 1\n")

# re_aa2dna function
#
#   RegEx, Dictionary -> RegEx
#   aaSeq  CodonTable    ntSeq
# this function translates a given regular expression amino acid sequence into
# a regular expression DNA sequence, given a specific genetic code table
# tests
# re_aa2dna('G?[AY]{2}.') --> (GG.)?[(GC.)(TA[TC])]{2}.{3}
# re_aa2dna('A(CE)+[^T]') -->
# re_aa2dna(G{1,3}.[^P]{2}[RGI]*L+V) -->
# --> (GG.){1,3}.{3}[^(CC.)]{2}[(CG.)(GG.)(AG[GA])(AT[TAC])]*[(CT.)(TT[GA])+(GT.)
#
# still doesn't support:
#  - any aa: '.' in or statement
#  - sequences of aa as patterns, eg. 'ARY' TODO
TASK_PATTERN = 'G{1,3}.[^P]{2}[RGI]*L+V'

def re_aa2dna(aaSequence = TASK_PATTERN,
              geneticCode = CodonTable.standard_dna_table.forward_table):


    #  List -> RegEx
    # This function accepts a list of codons, and reduces them to a regular
    # expression, where [ATCG] is . and three possible nucleotides are written
    # as the opposite of the remaining nucleotide [^X]
    # it only reduces codons in regard og their last nucleotide, taking from
    # the fact that it is usually the most variable one
    def minimizeSequence(listOfCodons):
        ntList = ['A','T','G','C']
        # print(listOfCodons)
        for codon in listOfCodons: codon = codon.upper()

        # create a list of generated regular expressions of nucleotides
        listOfCodonRegExp = []

        for nt1 in ntList:
            for nt2 in ntList:
                codonSubList = []
                # get all the codons that start with two specific nucleotides
                for codon in listOfCodons:
                    if codon[0] == nt1 and codon[1] == nt2:
                        codonSubList.append(codon)

                regexCodons = len(codonSubList)
                # depending on the number of codons starting with the same 2 nt
                # append list of codon regexp with a specific regexp
                if regexCodons == 0:
                    pass
                elif regexCodons == 1:
                    # for 1 codon just append it to the list
                    listOfCodonRegExp.append(codonSubList[0])
                elif regexCodons == 2:
                    # for 2 codons append them in form (nt1 nt2 [nt3_1 nt3_2])
                    nt3 = ''.join((codonSubList[0][2],codonSubList[1][2]))
                    listOfCodonRegExp.append(''.join((nt1,nt2,'[',nt3,']')))
                elif regexCodons == 3:
                    # for 3 codons append them in form (nt1 nt2 [^ not_nt3])
                    nt3 = [codonSubList[0][2],codonSubList[1][2],codonSubList[2][2]]
                    for nt in ntList:
                        if nt not in nt3:
                            not_nt3 = nt
                    listOfCodonRegExp.append(''.join((nt1,nt2,'[^',not_nt3,']')))
                elif regexCodons == 4:
                    # for 4 codons append them in form (nt1 nt2 .)
                    listOfCodonRegExp.append(''.join([nt1,nt2,'.']))
                else:
                    print('codon sub list is longer than 4, ')

        # print(listOfCodonRegExp)

        # create the final regular expression, by appending singular regexps
        codonRegExp = ""

        for regExp in listOfCodonRegExp:
            if len(listOfCodonRegExp) > 1:
                codonRegExp += regExp + '|'
            else:
                codonRegExp += '(' + regExp + ')'


        if len(listOfCodonRegExp) > 1:
            codonRegExp = '(' + codonRegExp[:len(codonRegExp)-1] + ')'

        # print(codonRegExp)
        return codonRegExp


    reverseGeneticCode = {}
    # reverse the genetic code dictionary
    for codon in geneticCode:
        if geneticCode[codon] not in reverseGeneticCode:
            reverseGeneticCode[geneticCode[codon]] = [codon]
        else:
            reverseGeneticCode[geneticCode[codon]].append(codon)

    # print(reverseGeneticCode)

    # regexp for finding aminoacids in the pattern
    aaPattern = '[A-Z]+'

    # find all aminoacid literals and list them in a non-redundant way
    listOfaa = set(re.findall(aaPattern, aaSequence))

    # a dictionary with aa regexps as keys, and a nt regexps as values
    aa2dnaDict = {}

    for literal in listOfaa:
        listOfCodons = []
        for i in range(len(literal)):
            listOfCodons.extend(reverseGeneticCode[literal[i]])
        aa2dnaDict[literal] = minimizeSequence(listOfCodons)

    ntSequence = aaSequence
    # first change all . occurrences to (...)
    ntSequence = re.sub('\.','(...)',ntSequence)

    for literal in listOfaa:
        # exchange all occurences of aminoacid regexps to nucleotide regexps
        # given the created dictionary aa2dnaDict
        ntSequence = re.sub('(?P<first>[^A-Z|\.]|^)' + literal +'(?P<last>[^A-Z|\.]|$)',
                            '\g<first>'+aa2dnaDict[literal]+'\g<last>',ntSequence)

    return ntSequence
'''
a) czy sekwencja aminokwasowa opisana danym wyrażeniem znajduje się
w sekwencji aminokwasowej hemoglobiny (uzyskanej w zadaniu 1)
'''
print("podpunkt a)")
print("Poszukiwanie sekwencji aminokwasowej odpowiadającej wyrażeniu regularnemu: %s \n" % TASK_PATTERN)

# load protein fasta file
try:
    hemoglobinProtein = SeqIO.read("hemoglobin_protein.fasta","fasta")
except FileNotFoundError:
    proteinFileName = input("Podaj ścieżkę do pliku fasta z sekwencją białka: ")
    while True:
        try:
            hemoglobinProtein = SeqIO.read(proteinFileName,"fasta")
            break
        except FileNotFoundError:
            mRNA_fileName = input("Coś jest nie tak, spróbuj jeszcze raz (ścieżka do pliku z mRNA): ")

#print(hemoglobinProtein, dir(hemoglobinProtein))
#print('\n',str(hemoglobinProtein.seq))
#print(TASK_PATTERN)
print('Znalezione sekwencje w białku:')
print(re.findall(TASK_PATTERN,str(hemoglobinProtein.seq)))

'''
b) czy sekwencja odpowiadająca równoważnemu wyrażeniu w języku nukleotydów
zostanie znaleziona w sekwencji nukleotydowej mRNA hemoglobiny
'''
print("\npodpunkt b)")
print("Poszukiwanie sekwencji nukleotydowej opisanej wyrażeniem regularnym "
      "odpowiadającym aminokwasowemu wyrażeniu regularnemu.")

# ntPattern = re_aa2dna(TASK_PATTERN)
# print(ntPattern)
# manually editing pattern, function didn't succeed
ntPattern = "(GG.){1,3}...([AGT]..|.[AGT].){2}(AT^G|AGA|AGG|GG.|CG.)*(TTG|TTA|CT.)+GT." # courtesy of Jacek Śmietański
print(ntPattern,'\n')

# load mRNA fasta file
try:
    hemoglobin_mRNA = SeqIO.read("hemoglobin_mRNA.fasta","fasta")
except FileNotFoundError:
    mRNA_fileName = input("Podaj ścieżkę do pliku fasta z sekwencją mRNA: ")
    while True:
        try:
            hemoglobin_mRNA = SeqIO.read(mRNA_fileName,"fasta")
            break
        except FileNotFoundError:
            mRNA_fileName = input("Coś jest nie tak, spróbuj jeszcze raz (ścieżka do pliku z mRNA): ")

hemoglobin_mRNA = SeqIO.read("/Users/Kozel/Documents/UJ/Biotechnologia molekuarna/3 semestr/Bioinformatyka/Kody/hemoglobin_mRNA.fasta","fasta")

for match in re.finditer(ntPattern,str(hemoglobin_mRNA.seq)):
    s = match.start()
    e = match.end()
    print("W mRNA znaleziono dopasowanie: %s na pozycjach od %d do %d" %(hemoglobin_mRNA.seq[s:e], s, e))
print("\nJest to sekwencja, która koduje odnaleziony wcześniej fragment białka.")

'''
część 2:
Pewna charakterystyczna grupa białek, która posiada zdolność wiązania się
do helisy DNA ma w swojej sekwencji następujący motyw:

([^.{0,3}]| Q).[^FHWY][ILM][^P][^FHILVWYP][DHFM][FMY]

Sprawdź, które z sekwencji z pliku amino_acid_sequences.txt mogą być
fragmentami łańcucha z tej rodziny białek oraz określ, który z trzech
aminokwasów I, L czy M pojawiał się najczęściej w miejscu sekwencji
określonym wyrażeniem [ILM]
'''
#
print("\n\nCzęść 2\n")

print('Poszukiwanie fragmentów białek z danej rodziny.\n')
TASK2_PATTERN = "([^.{0,3}]| Q).[^FHWY](?P<query>[ILM])[^P][^FHILVWYP][DHFM][FMY]"

try:
    f = open("amino_acid_sequences.txt",'r')
except FileNotFoundError:
    aaFileName = input("Podaj ścieżkę do pliku z sekwencjami aminokwasów: ")
    while True:
        try:
            f = open(aaFileName,'r')
            break
        except FileNotFoundError:
            aaFileName = input("Podaj ścieżkę do pliku z sekwencjami aminokwasów raz jeszcze: ")

aa_sequences = f.read()
f.close()

whichAminoAcid = {'I':0,'L':0,'M':0}

for match in re.finditer(TASK2_PATTERN,aa_sequences):
    s = match.start()
    e = match.end()
    line = aa_sequences[0:e].count('\n')
    whichAminoAcid[match.group('query')] += 1
    print("Znaleziono %s na pozycjach od %d do %d, jest to sekwencja %d." %(hemoglobin_mRNA.seq[s:e], s, e, line+1))

maxaaOccurences = max(whichAminoAcid.values())
mostFrequentaa = []
for aa in whichAminoAcid:
    if whichAminoAcid[aa] == maxaaOccurences:
        mostFrequentaa.append(aa)

print('\nNajczęściej występującym aminokwasem w pozycji określonej wyrażeniem '
      '[ILM] jest %s.' % ' i '.join(mostFrequentaa))
