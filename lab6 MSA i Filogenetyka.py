__author__ = 'Kozel'
# jak chce się robić bardziej uniwersalne programy, to się korzysta z tego
# modułu do wczytywania plików etc.
# import os.path
from Bio import AlignIO
from Bio import Phylo


alignment = AlignIO.read("rho.clustalw","clustal")
print(alignment)

# wypisanie w formacie fasta
print("\nFormat fasta")
print(alignment.format("fasta"))

# sekwencja i id
for record in alignment:
    print(record.seq, record.id)

# wypisanie zadanych elementów, czyli sekwencji
print(alignment[4:7])
# same wybrane residues dla wybranej sekwencji
print(alignment[3].seq[6:20])
# konkretne residues dla wszystkich sekwencji
print(alignment[:,1:6])

###### drzewa filogenetyczne

print("\n\nDrzewa filogenetyczne\n")
tree = Phylo.read("rhoTree.ph", "newick")
print(tree)

# wizualizacja znakami ASCII
Phylo.draw_ascii(tree)

# Phylo.draw_graphviz(tree)

# operacje na drzewie
terminals = tree.get_terminals()
for terminal in terminals:
    print(terminal)

# droga do człeka od korzenia(?)
print("\nPath do człeka:")
path = tree.get_path("Człek")
print(path)

print("\nWspólny przodek człowieka i dania:")
common_ancestor = tree.common_ancestor("Człek", "Danio")
print(common_ancestor)

# głebokości w drzewie
print("Depths in the tree:")
depths = tree.depths()
# print(depths)
for d in depths.keys():
    print(d, ":", depths.get(d))
