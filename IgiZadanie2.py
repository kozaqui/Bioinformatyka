__author__ = 'Iga'
from Bio import Entrez
Entrez.email = "iga.niemiec@student.uj.edu.pl"

handle = Entrez.einfo(db="nucleotide")
result = Entrez.read(handle)
print(result["DbInfo"].keys())

for field in result["DbInfo"] ["FieldList"]:
    print("%(Name)s, %(FullName)s, %(Description)s" % field)


handle = Entrez.esearch(db="nucleotide", term="BRCA1[TITL] AND Homo sapiens[ORGN]", retmax="100")
result = Entrez.read(handle)
list1 = result["IdList"]
print(list1)

handle2 = Entrez.esummary(db="nucleotide", id=",".join(list1))
result2 = Entrez.read(handle2)

for record in result2:
    print ("%(Name)s, %(Chromosome)s, %(MapLocation)s, %(Description)s" % record)
'''
for record in list1:
    handle = Entrez.efetch(db="gene", id=record, rettype="gb", retmode="text")
    print(handle.read())
'''
###################################
# BRCA1 ID=672
'''
handle = Entrez.einfo(db="omim")
result = Entrez.read(handle)
print(result["DbInfo"].keys())

for field in result["DbInfo"] ["FieldList"]:
    print("%(Name)s, %(FullName)s, %(Description)s" % field)

handle2 = Entrez.esearch(db="omim", term="672[ID]", retmax="100")
result = Entrez.read(handle2)
list2 = result["IdList"]
print(list2)

'''