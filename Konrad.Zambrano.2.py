from Bio import Entrez
from Bio import SeqIO

Entrez.email = "konrad.zambrano_quispe@student.uj.edu.pl"

def ELinkToList(linkHandle):
# this function accepts a handle to
    resultLink = Entrez.read(linkHandle)

    linkList = []
    for link in resultLink[0]['LinkSetDb'][0]['Link']:
        linkList.append(link["Id"])

    return linkList

'''
Napisz
skrypt, który przy wykonuje kolejno następujące zadania:
1. Wyszukaj geny powiązane z nazwą BRCA1 u człowieka. Wypisz pola: Name,
   Chromosome, MapLocation i Description
'''
print('\nZadanie 1\n')

handle = Entrez.esearch(db="gene", term = "BRCA1 AND Homo sapiens[ORGN]", retmax = 10)
resultSearch = Entrez.read(handle)
# print(resultSearch)
# {'RetMax': '20', 'TranslationSet': [], 'RetStart': '0', 'Count': '12419', 'IdList': ['7157', '1956', '7422', '4524', '22059', '7040', '2099', '155971', '672', '2064', '3091', '4790', '367', '5243', '1029', '207', '1312', '7421', '3845', '6774'], 'TranslationStack': [{'Field': 'All Fields', 'Explode': 'N', 'Term': 'BRCA1[All Fields]', 'Count': '12419'}, 'GROUP'], 'QueryTranslation': 'BRCA1[All Fields]'}
handle = Entrez.esummary(db="gene", id=",".join(resultSearch["IdList"]))
result = Entrez.read(handle)

print('Ludzkie geny powiązane z nazwą BRCA1:')
i = 0
for record in result['DocumentSummarySet']['DocumentSummary']:
    print("%(Name)s, %(Chromosome)s, %(MapLocation)s, %(Description)s" % record)
    if record['Name'] == "BRCA1":
        geneId = resultSearch["IdList"][i]
    i += 1
#result = Entrez.read(handle)
'''
2. Zachowaj id genu o nazwie BRCA1
'''
print('\n\nZadanie 2\n')

print('Id genu BRCA1 to:')
print(geneId)
# geneId = 672
'''
3. Wypisz tytuły rekordów bazy OMIM związane z tym genem
'''
print('\n\nZadanie 3\n')

handleLink = Entrez.elink(dbfrom="gene", db="omim", id = geneId)

omimLinks = ELinkToList(handleLink)


handleOmim = Entrez.esummary(db="omim", id=",".join(omimLinks))
resultOmim = Entrez.read(handleOmim)

print()

print("Tytuły rekordów z bazy OMIM, powiązane z genem BRCA1:")
for record in resultOmim:
    print(record['Title'])

'''
4. Znajdź białka powiązane z tym genem (baza protein): wypisz pola numeru
   identyfikacyjnego Gi dla wyszukanych sekwencji
'''
print('\n\nZadanie 4\n')

handleLink = Entrez.elink(dbfrom="gene", db="protein", id = geneId)
proteinLinks = ELinkToList(handleLink)

handleProtein = Entrez.esummary(db="protein", id=",".join(proteinLinks))
resultProtein = Entrez.read(handleProtein)

# print(resultProtein)
print()

print("Identyfikatory Gi dla białek powiązanych z genem BRCA1")
for record in resultProtein:
    print("%(Gi)s" % record)
'''
5. Pobierz białko o numerze identyfikacyjnym 121949022 w formacie genbank,
   wypisz jego nazwę oraz sekwencję aminokwasów
'''
print('\n\nZadanie 5\n')

handle = Entrez.efetch(db="protein", id='121949022', rettype="gb",
                       retmode="text")
result = SeqIO.read(handle, "genbank")
print(handle.read())
handle.close()

print('Nazwa rekordu białka o id: 121949022')
print(result.name)
print('Sekwencja białka:')
print(result.seq)
print()
'''
6. Znajdź wszystkie mutacje w obrębie genu BRCA1 u człowieka
   (w elink do bazy SNP użyj term="Homo sapiens"), wypisz ile
   mutacji znaleziono, dla pierwszych 50 wyświetl ich numer
   identyfikacyjny SNP_ID, klasę SNP_CLASS, gen GENE i pozycję CONTIGPOS
'''
print('\n\nZadanie 6\n')

handleLink = Entrez.elink(dbfrom='gene', db='snp', id=geneId,
                          term="Homo sapiens")
snpList = ELinkToList(handleLink)

print()
print("Odnaleziono %d mutacji dla genu BRCA1" % len(snpList))

handle = Entrez.esummary(db='snp',id=",".join(snpList[:50]))
result = Entrez.read(handle)

print('50 pierwszych rekordów SNP genu BRCA1 u człowieka:')
for record in result:
    print("SNP ID: %(SNP_ID)s, klasa SNP: %(SNP_CLASS)s, gen: %(GENE)s"
          ", CONTIGPOS: %(CONTIGPOS)s." % record)
'''
7. Na podstawie przeprowadzonych wyżej poszukiwań, sformułuj odpowiednie
   wnioski (jakie informacje udało Ci się uzyskać? jakie jest ich znaczenie?).
   Treść wniosku wyświetl na ekranie (print) lub zawrzyj w komentarzu na końcu
   pliku ze skryptem.
'''
print('\n\nZadanie 7\n')

print('''
    Pierwszy wniosek: gen BRCA1 jest bardzo dokładnie przebadanym genem,
    jest o nim całe mnóstwo informacji w bazie NCBI.

    Zadanie 1
    Powiązane z nim są również geny czynników wzrostowych oraz wiele
    innych genów związanych z transformacją nowotworową.

    Zadanie 3
    Sugeruje, że mutacje w genie BRCA1 mogą się wiązać z rodzinnymi
    skłonnościami do zapadania na nowotwory piersi i jajnika, a także
    na podatność na raka trzustki.
    Nie dziwota, że Angelina Jolie zdecydowała się na mastektomię.

    Zadanie 4
    Jest bardzo dużo białek powiązanych z tym genem.

    Zadanie 6
    Być może wynika to z jeszcze większej liczby mutacji, jakie zostały dla tego
    genu opisane. Ilość tych mutacji może wskazywać na to, że ten region genomu
    może wykazywać większe tempo mutacji, lub ta duża ilość opisanych mutacji
    jest pochodną tego, że gen ten jest ekstensywnie badany.
    ''')