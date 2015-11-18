__author__ = 'Kozel'
# encoding utf-8
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.Emboss import Applications

##### wczytuję primery i egzony konstruktu
fragments = list(SeqIO.parse("/Users/Kozel/Documents/UJ/Golec project/Plan działania/2.2 H_sapiens_HAS2_construct_fragments.fasta","fasta"))
primers = list(SeqIO.parse("/Users/Kozel/Documents/UJ/Golec project/Plan działania/primery.fasta","fasta",alphabet=IUPAC.unambiguous_dna))

for primer in primers:
    print("Description: %s" % primer.description)
    print("Sequence: %s" % primer.seq)
    print("Complement: %s" % primer.seq.complement())
    print("Reverse complement: %s\n" % primer.seq.reverse_complement())

## teraz restriction
# from Bio import Restriction
# print(Restriction.AllEnzymes)
# print(Restriction.AgeI.site)

################################################
##### tworzę sekwencję całego wektora ze wstawką
################################################

vector = SeqIO.read("/Users/Kozel/Documents/UJ/Golec project/Wektory/pLX304.txt",
                    'fasta',IUPAC.unambiguous_dna)
insert = SeqIO.read("/Users/Kozel/Documents/UJ/Golec project/Wektory/pLX304_HAS2_insert.txt",
                    'fasta',IUPAC.unambiguous_dna)

print("\n\nInfo o wektorze:")
print(vector)
print(dir(vector))
print(vector.features)

print(dir(SeqFeature))

my_start_pos = SeqFeature.ExactPosition(2)
my_end_pos = SeqFeature.ExactPosition(24)
my_feature_location = SeqFeature.FeatureLocation(my_start_pos,my_end_pos)
my_feature_type = "CDS" # "CDS"

# my_feature = SeqFeature(my_feature_location,type=my_feature_type)
my_feature = SeqFeature(my_feature_location,type=my_feature_type)

vector.features.append(my_feature)

print(vector.features)

# Applications.Primer3Commandline