from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import re

sequence = Seq("ATGtcccacta",IUPAC.unambiguous_dna)
print(sequence)

#4 nić komplementarna
print(sequence.complement())
# ::-1 to to samo co [0:len(seq):-1], początek, koniec i krok, sam : oznacza domyślne [0:len(seq):-1]
print(sequence.complement()[::-1])
# to samo co
print(sequence.reverse_complement())

simple_seq = Seq("TCTGTGCTAAAGTGTAACTCGTAGGCACTATCTAC")
simple_seq_r = SeqRecord(simple_seq, id="AC834343", name="seqX", description="Homo erectus, chr25")
print(simple_seq_r)

record = SeqIO.read("/Users/Kozel/Documents/UJ/Biotechnologia molekuarna/3 semestr/Bioinformatyka/Kody/hemoglobin.txt","fasta")
print(record.id)
print(record.name)
print(record.description)
print(record.annotations)
#print(record.seq)

## wyrażenia regularne
#re.match(pattern,string) - szuka na początku stringa
#re.search(pattern,string) - szuka w całym stringu
print(re.search("C.T","ATCATGGC"))

# często daje się na początek ^ a na koniec $

## zadanie 7
print('\nZadanie 7\n')

seq_aa ="MVHLTPEEKSAVTALW"
pattern1 = "^[^P]+T*[AT].{3,5}..[VCPGH]C?.*W$"
match = re.search(pattern1, seq_aa)
print(match)
print(match.start())
print(match.end())

## zadanie 8
print('\nZadanie 8\n')

pattern2 = "^[^P]+T*[AT].{3,5}(..)[VCPGH]C?.*W$"
match = re.search(pattern2, seq_aa)
print(match.group(0))
print(match.group(1))

#
# for i in range(len(match.group())):
#     print(i)
#     print(match.group(i))

## zadanie 9
print('\nZadanie 9\n')

# nie tworzy grupy
pattern3 = "M[^T]H(?=LTP)(.*)"
match = re.search(pattern3, seq_aa)
print(match.group(0))
print(match.group(1))

pattern4 = "H(?=LTV)" # nie będzie dopasowania
match = re.search(pattern4, seq_aa)
print(match)

## trying

print('\n\ntrying:')
string = "TYPOTTTYUPOY"
#re.findall(r'(?=(\w\w))', 'hello')
#['he', 'el', 'll', 'lo']
pattern = "(?=(T|P))(T|P)*"
#print(re.search(pattern,string))
print(re.findall(pattern,string))

for match in re.finditer(pattern,string):
    s = match.start()
    e = match.end()
    print("Znaleziono %s na pozycjach od %d do %d." %(string[s:e], s, e,))
# import this
# this.c