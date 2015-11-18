__author__ = 'Kozel'
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

print(help(NCBIWWW.qblast))

# result_handle = NCBIWWW.qblast("blastn", "nr", "23527284")
# print(result_handle.read())

# s = SeqIO.read(open("sequence.fa"), format="fasta")
# result_handle = NCBIWWW.qblast("blastn", "nt", s.seq)

# zmniejszanie liczby wyników za pomocą hitlist_size
s = SeqIO.read("lab5-seq1.fasta", format="fasta")
result_handle = NCBIWWW.qblast("blastn", "nx", s.seq, hitlist_size=1)

# bazy:
# nr - sekwencje białkowe
# nt - nukleotydowe
# inne:
# blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE= BlastDocs&DOC_TYPE =ProgSelectionGuide#db$)$

blast_records = NCBIXML.parse(result_handle)
blast_record = blast_records.next()

for alignment in blast_record.alignments:
    print('Alignment------------------------------')
    print('title: ', alignment.title)
    print('length: ', alignment.length)
    for hsp in alignment.hsps:
        print('HSP : ')
        print('e value:', hsp.excpect)
        print(hsp.query[0:75] + '...')
        print(hsp.match[0:75] + '...')
        print(hsp.sbjct[0:75] + '...')