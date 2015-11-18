__author__ = 'Kozel'
# dictionary.items()

from Bio.Data import CodonTable
import re
import sre_parse

from sys import argv
#script, fileIn = argv

aa = 'R'
geneticCode = CodonTable.standard_dna_table.forward_table
#print(geneticCode)
for codon in geneticCode:
    if geneticCode[codon] == aa:
        print(codon)

def re_aa2dna(aaSequence, geneticCode = CodonTable.CodonTable.forward_table):
    return '[ATGC]'

# re_aa2dna('a')

#   RegEx -> RegEx
# this function translates a given regular expression amino acid sequence into a regular expression DNA sequence
# given a specific genetic code table
# testuj
# re_aa2dna('G?[AY]{2}.') --> (GG.)?[(GC.)(TA[TC])]{2}.{3}
# re_aa2dna('A(CE)+[^T]') -->
# re_aa2dna(G{1,3}.[^P]{2}[RGI]*L+V) --> (GG.){1,3}.{3}(?!=CC.){2}[(CG.)(GG.)(AG[GA])(AT[TAC])]*[(CT.)(TT[GA])+(GT.)
#

#def proteinToDNA(aaSequence, geneticCode = CodonTable.CodonTable.forward_table)
'''
gencode = { 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            ‘CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W' }
'''

''' Potrzebne są:
    funkcja, która szuka jednego elementu regexp
    funkcja, która zamienia dany element z aa na kodon
    funkcja, która grupuje kodony w dodatkowy regexp

    Opcjonalnie:
    funkcja, która sprawdza poprawność inputu
     . . . . . .
    Inny sposób implementacji:
    funkcja, która szuka re ciągu aminokwasów, a także dowolnego aminokwasu, czyli kropki
    funkcja, która wypisuje wszystkie możliwe kodony i redukuje je, jeśli jest taka możliwość
    funkcja która zamienia wejsciowe elementy na obliczone el. wyjsciowe
'''
patterns = ['G?[AY]{2}.', 'A(CE)+[^T]', 'A*', '.','.*','ACD']
test_pattern = 'A(CE)+[^T]'

#re.match()
import _sre
import sre_constants
import sre_parse
import sre_compile
from sre_parse import Tokenizer, Pattern, parse_template, parse
print('\n\n')

#print(parse(pattern))
# for element in parse(pattern):
#     print(element)

for pattern in patterns:
    print('\n')
    print(pattern)
    for element in parse(pattern):
        for subelement in element:
            print(subelement)
        print(element)
#parse_template('asdadasda', pattern)

# print(dir(Tokenizer(test_pattern)))
# k = Tokenizer(test_pattern)
# print(k.index)
# print(k.get)
# k.index = 3
# print(k.index)
# print(k.next)
# print(Tokenizer(test_pattern).next)
# print(Tokenizer(test_pattern).next)

# pattern = "[abc][a-z]*.?[AHF]{2,5}"
# prog = re.compile(pattern)
# print('A oto pattern: %s i skompilowany pattern: %s.' % (pattern, prog.pattern))
# newThing = sre_parse.parse(pattern)
# print(newThing)
# print("")
# print(newThing.data)
print(chr(65))

reverseGeneticCode = {}
# odwróć słownik
for codon in geneticCode:
    if geneticCode[codon] in reverseGeneticCode:
        reverseGeneticCode[geneticCode[codon]].append(codon)
    else:
        reverseGeneticCode[geneticCode[codon]] = [codon]

print(reverseGeneticCode)
print(reverseGeneticCode.keys())
# print(geneticCode.values())
# for aa in geneticCode.values():
#     for
#     reverseGeneticCode[aa] =
# szuknie patternów aminokwasów
aaPattern = '[A-Z|\.]+' # znajdź w regexp symbol aminokwasu/ów
taskPattern = 'G{1,3}.[^P]{2}[RGI]*L+V'


print(re.findall(aaPattern, taskPattern))
print(re.finditer(aaPattern,taskPattern))

# znajduje wszystkie regex aminokwasy
for match in re.finditer(aaPattern,taskPattern):
    s = match.start()
    e = match.end()
    print("Znaleziono %s na pozycjach od %d do %d" %(taskPattern[s:e], s, e))


# print(re.match(aaPattern, taskPattern))
# print(re.search(aaPattern, taskPattern))

# regexp do zastępowania patterna innym
#text = re.sub()
aa2nt = re.compile(aaPattern)
ntSequence = aa2nt.sub('(ATG)',taskPattern)
print(ntSequence)

#


# bold = re.compile(r'\*{2}(.*?)\*{2}', re.UNICODE)
# text = 'Make this **bold**. This **too**.'
#
# print('Text:', text)
# print('Bold:', bold.sub(r'<b>\1</b>', text))# , count=1))