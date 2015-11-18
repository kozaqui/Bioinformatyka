__author__ = 'Kozel'
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "konrad.zambrano_quipe@student.uj.edu.pl" # Dżem dobry

## zadanie 2

## po przedstawieniu się
# handle = Entrez.einfo()
# result = handle.read()
# print("Wersja w formacie XML:")
# print(result)
# For some reason it has to be loaded for a second time
handle = Entrez.einfo()
result = Entrez.read(handle) # używając parsera z Bio.Entrez
print("Wersja korzystająca z parsera Bio.Entrez:")
print(result)

## zadanie 3

# zobacz informacje o konkretnej bazie (Nucleotide:
# handle = Entrez.einfo(db='nucleotide')
# result = handle.read()
# print("\n\nBaza Nucleotide. Wersja w formacie XML:")
# print(result)

# handle = Entrez.einfo(db='nucleotide')
# result = Entrez.read(handle)
# print("\n\nInfo o bazie Nucleotide:")
# print(result)
#
# print("\nKeys")
# print(result["DbInfo"].keys())
# print("\nDescription")
# print(result["DbInfo"]["Description"])
# print("\nCount")
# print(result["DbInfo"]["Count"])
# print("\nLastUpdate")
# print(result["DbInfo"]["LastUpdate"])
#
# print("Field List")
# print(result["DbInfo"]["FieldList"])
# for field in result["DbInfo"]["FieldList"]: # pola, po których możemy wyszukiwać
#     print("%(Name)s, %(FullName)s, %(Description)s" % field)

# Sprawdzanie innych baz
# handle = Entrez.einfo(db='nucgss')
# handle = Entrez.einfo(db='popset') # sekwencje do badania drzew filogenetycznych

## Zadanie 4

# handle = Entrez.esearch(db="nucleotide", term="swine influenza", retmax = "20")
# result = Entrez.read(handle)
# print("\nResult ID of nucleotides connected to  list of length: %d" % len(result["IdList"]))
# print(result["IdList"])
#
term = "influenza AND H1N1[TITL] AND (hemagglutinin[PROT] OR neuraminidase[PROT])"
handle = Entrez.esearch(db="nucleotide", term=term, retmax = "20")
result = Entrez.read(handle)
# print("\nResult ID of nucleotides connected to  list of length: %d" % len(result["IdList"]))
# print(result["IdList"])

# handle = Entrez.esearch(db="popset", term = term, retmax = "20")
# result = Entrez.read(handle)
# print(result)

## Zadanie 5
# handle = Entrez.esummary(db="nucleotide", id="326579398") # pojedynczy rekord
#
# print(result["IdList"])
# handle = Entrez.esummary(db="nucleotide", id=",".join(result["IdList"]))
# result = Entrez.read(handle)
#
# for record in result:
#     print("%(Title)s, %(Gi)s, %(Length)s %(CreateDate)s" % record)
#
# ## Zadanie 6
#
# # pobieranie w innym formacie, tutaj genbank
# handle = Entrez.efetch(db="nucleotide", id="326579398", rettype = "gb", retmode = "text")
# print(handle.read())
#
# # wczytanie do obiektu SeqRecord
# handle = Entrez.efetch(db="nucleotide", id="326579398", rettype="gb", retmode = "text")
# record = SeqIO.read(handle, "genbank")
# print(record.name)
# print(record.description)
# print(record.seq)

## Zadanie 7
# handle = Entrez.esearch(db="protein", term="Homo sapiens[orgn] and rhodopsin",
#                         usehistory = "y")
# result = Entrez.read(handle)
# query_key = result['QueryKey']
# webenv = result['WebEnv']
#
# handle = Entrez.efetch(db="protein", WebEnv = webenv, query_key = query_key,
#                        retstart = 0, retmax = 2, rettype = "gp", retmode = "text")
# print(handle.read())

## Zadanie 8

# ID genu hemoglobiny
handle = Entrez.elink(dbfrom = "gene", db="protein", id=3043)
result = Entrez.read(handle)

print(result)

# dla pierwszego id i pierwszej bazy docelowej
for link in result[0]["LinkSetDb"][0]["Link"] :
    print(link["Id"]) # lista linków

# form Bio.Blast.Applications