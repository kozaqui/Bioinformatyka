__author__ = 'Kozel'
from Bio.Data import CodonTable
import re

# re_aa2dna('a')

#   RegEx -> RegEx
# this function translates a given regular expression amino acid sequence into
# a regular expression DNA sequence, given a specific genetic code table
# tests
# re_aa2dna('G?[AY]{2}.') --> (GG.)?[(GC.)(TA[TC])]{2}.{3}
# re_aa2dna('A(CE)+[^T]') -->
# re_aa2dna(G{1,3}.[^P]{2}[RGI]*L+V) --> (GG.){1,3}.{3}[^(CC.)]{2}[(CG.)(GG.)(AG[GA])(AT[TAC])]*[(CT.)(TT[GA])+(GT.)
#
# still doesn't support:
#  - any aa: '.' in or statement
#  - sequences of aa as patterns, eg. 'ARY' TODO
TASK_PATTERN = 'G{1,3}.[^P]{2}[RGI]*L+V'

#   RegEx -> RegEx
# this function translates a given regular expression amino acid sequence into
# a regular expression DNA sequence, given a specific genetic code table
# tests
# re_aa2dna('G?[AY]{2}.') --> (GG.)?[(GC.)(TA[TC])]{2}.{3}
# re_aa2dna('A(CE)+[^T]') -->
# re_aa2dna(G{1,3}.[^P]{2}[RGI]*L+V) --> (GG.){1,3}.{3}[^(CC.)]{2}[(CG.)(GG.)(AG[GA])(AT[TAC])]*[(CT.)(TT[GA])+(GT.)
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

        listOfCodonRegExp = []

        for nt1 in ntList:
            for nt2 in ntList:
                codonSubList = []
                # get all the codons that start with two specific nucleotide
                for codon in listOfCodons:
                    if codon[0] == nt1 and codon[1] == nt2:
                        codonSubList.append(codon)

                regexCodons = len(codonSubList)
                if regexCodons == 0:
                    pass
                elif regexCodons == 1:
                    listOfCodonRegExp.append(codonSubList[0])
                elif regexCodons == 2:
                    nt3 = ''.join((codonSubList[0][2],codonSubList[1][2]))
                    listOfCodonRegExp.append(''.join((nt1,nt2,'[',nt3,']')))
                elif regexCodons == 3:
                    nt3 = [codonSubList[0][2],codonSubList[1][2],codonSubList[2][2]]
                    for nt in ntList:
                        if nt not in nt3:
                            notnt3 = nt
                    listOfCodonRegExp.append(''.join((nt1,nt2,'[^',notnt3,']')))
                elif regexCodons == 4:
                    listOfCodonRegExp.append(''.join([nt1,nt2,'.']))
                else:
                    print('codon sub list is longer than 4, ')

        # print(listOfCodonRegExp)
        codonRegExp = ""

        for regExp in listOfCodonRegExp:
            codonRegExp += '(' + regExp + ')'

        if len(listOfCodonRegExp) > 1:
            codonRegExp = '[' + codonRegExp + ']'

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
    # aaPattern = '[A-Z]+'

    # find all aminoacid literals and list them in a non-redundant way
    listOfaa = set(re.findall(aaPattern, aaSequence))
    # print(listOfaa)

    aa2dnaDict = {'\.':'(...)'}

    for literal in listOfaa:
        listOfCodons = []
        for i in range(len(literal)):
            listOfCodons.extend(reverseGeneticCode[literal[i]])
        # print('\n')
        # print(literal)
        aa2dnaDict[literal] = minimizeSequence(listOfCodons)

    ntSequence = aaSequence
    ntSequence = re.sub('\.','(...)',ntSequence)
    #aa2nt = re.compile(aaPattern)
    for literal in listOfaa:
        #print(literal)
        #print(aa2dnaDict[literal])
        ntSequence = re.sub('(?P<first>[^A-Z|\.]|^)' + literal +'(?P<last>[^A-Z|\.]|$)','\g<first>'+aa2dnaDict[literal]+'\g<last>',ntSequence)
        #print(ntSequence)

# for pattern in [ r’^(?P<first_word>\w+)’, r’(?P<last_word>\w+)\S*$’,
# r’(?P<t_word>\bt\w+)\W+(?P<other_word>\w+)’, r’(?P<ends_with_t>\w+t)\b’,
# ]:


    #print(aaSequence)
    #print(ntSequence)

    return ntSequence


re_aa2dna()
print(re_aa2dna('G?[AY]{2}.'))
print(re_aa2dna('A(CE)+[^T]'))

''' Potrzebne są:
    funkcja, która szuka jednego elementu regexp
    funkcja, która zamienia dany element z aa na kodon
    funkcja, która grupuje kodony w dodatkowy regexp

    Opcjonalnie:
    funkcja, która sprawdza poprawność inputu
     . . . . . .
    Inny sposób implementacji:
    funkcja, która szuka re ciągu aminokwasów, a także dowolnego aminokwasu, czyli kropki - ok
    funkcja, która wypisuje wszystkie możliwe kodony - ok
    funkcja, która redukuje możliwe kodony do regexp
    funkcja która zamienia wejsciowe elementy na obliczone el. wyjsciowe - ok
'''

# reverseGeneticCode = {}
# # odwróć słownik
# for codon in geneticCode:
#     if geneticCode[codon] in reverseGeneticCode:
#         reverseGeneticCode[geneticCode[codon]].append(codon)
#     else:
#         reverseGeneticCode[geneticCode[codon]] = [codon]
#
# print(re.findall(aaPattern, taskPattern))
#
# for match in re.finditer(aaPattern,taskPattern):
#     s = match.start()
#     e = match.end()
#     print("Znaleziono %s na pozycjach od %d do %d" %(taskPattern[s:e], s, e))
#
# aa2nt = re.compile(aaPattern)
# ntSequence = aa2nt.sub('(ATG)',taskPattern)
# print(ntSequence)