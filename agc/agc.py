#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw
import operator

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

""" Dé-duplication en séquence “complète”  """

# Prend deux arguments correspondant au fichier fasta.gz et à la longueur minimale des séquences
# et retourne un générateur de séquences de longueur l >= minseqlen

def read_fasta(amplicon_file, minseqlen):

    lines = [l.rstrip('\n') for l in gzip.open(amplicon_file, 'rt')]
    sequence = ""
    for line in lines:
        if ">" not in line:
            sequence += line
        else:
            if len(sequence) >= minseqlen:
                yield sequence
            sequence = ""

# Retourne un générateur des séquences uniques ayant une occurrence O>=mincount ainsi que leur occurrence

def dereplication_fulllength(amplicon_file, minseqlen, mincount):

    sequences_list = [seq for seq in read_fasta(amplicon_file, minseqlen)]
    dict = {seq_key:sequences_list.count(seq_key) for seq_key in sequences_list}
    sorted_dict = sorted(dict.items(), key=operator.itemgetter(1), reverse=True)
    for s in range(0, len(sorted_dict)):
        count = sorted_dict[s][1]
        if count >= mincount:
            yield [sorted_dict[s][0], count]

""" Recherche de séquences chimériques par approche “de novo” """

# Prend une séquence et une longueur de segment l: chunk_size et retourne une liste de sous-séquences de
# taille l non chevauchantes. A minima 4 segments doivent être obtenus par séquence

def get_chunks(sequence, chunk_size):

    sub_seq_list = []
    i = 0 #compteur
    s = chunk_size #variable qui conditionne la continuité de la boucle
    while(s <= len(sequence)):
        sub_sequence = sequence[i:i+chunk_size]
        sub_seq_list.append(sub_sequence)
        i += chunk_size
        s = i + chunk_size
    if len(sub_seq_list) >= 4:
        return sub_seq_list
    else:
        raise ValueError(" La taille de la séquence ne permet pas d'obtenir plus de 3 sous-séquences")

def get_unique(ids):
    return {}.fromkeys(ids).keys()

def common(lst1, lst2):
    return list(set(lst1) & set(lst2))

# Prend une séquence et une longueur de k-mer et retourne un générateur de tous les mots de longueur k présents dans cette séquence

def cut_kmer(sequence, kmer_size):

    length = len(sequence)
    for i in range(length-kmer_size+1):
        yield sequence[i:i+kmer_size]

# Prend un alignement (sous forme de liste) et calcule le pourcentage d’identité entre les deux séquences

def get_identity(alignment_list):

    seq1, seq2 = (alignment_list[0], alignment_list[1])
    id = sum(1 if c1 == c2 else 0 for c1, c2 in zip(seq1, seq2))
    return id / len(alignment_list[0]) * 100

# Prend un dictionnaire ayant pour clé un dictionnaire de kmer (vide au départ) et pour valeur une liste d’identifiant des séquences
# dont ils proviennent. Et retourne un dictionnaire de kmer contenant les kmers uniques présents dans chaque séquence pour une longueur de donnée de kmer

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):

    list_kmers = [kmer for kmer in cut_kmer(sequence, kmer_size)]
    for kmer in list_kmers:
        if kmer in kmer_dict.keys():
            kmer_dict[kmer].append(id_seq)
        else:
            kmer_dict[kmer] = [id_seq]
    return kmer_dict

# Prend un dictionnaire ayant pour clé un index de kmer et pour valeur une liste d’identifiant des séquences dont ils proviennent,
# une séquence et une longueur de kmer

def search_mates(kmer_dict, sequence, kmer_size):

    list_kmers = [kmer for kmer in cut_kmer(sequence, kmer_size)]
    c = Counter()
    for kmer in kmer_dict.keys():
        if kmer in list_kmers:
            c += Counter(kmer_dict[kmer])
    list_similar_sqs = list(list(zip(*(c.most_common(8))))[0])
    return list_similar_sqs

# prend une matrice donnant par segment le taux d’identité entre la séquence candidate et deux séquences parentes et retourne un booléen
# indiquant si la séquence candidate est une chimère (True) ou ne l’est pas (False)

def detect_chimera(perc_identity_matrix):

    seq0_similarity = False
    seq1_similarity = False
    std = []
    for elt in perc_identity_matrix:
        std.append(statistics.stdev(elt))
        if elt[0] > elt[1]:
            seq0_similarity = True
        elif elt[0] < elt[1]:
            seq1_similarity = True
    std = statistics.mean(std)
    return bool(std > 5 and seq0_similarity and seq1_similarity)

# Fait appel au générateur fourni par dereplication_fulllength et retourne un générateur
# des séquences non chimérique au format: yield [sequence, count]

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):

    sqs_count = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    sqs = list(map(itemgetter(0), sqs_count))
    perc_identity_matrix = []
    non_chimeres = list()
    id_seq = 0
    kmer_dict = {}
    for sequence, count in sqs_count:
        chunks = []
        # diviser chaque séquence candidate en 4 segments de longueur L= chunk_size
        segments = get_chunks(sequence, chunk_size)
        for seg in segments:
            mates_sqs = search_mates(kmer_dict, seg, kmer_size)
            chunks.append(mates_sqs)
        c = []
        for i in range(len(chunks)):
            c = common(c, chunks[i])
        if len(com) >= 2:
            for elt in c[0:2]:
                s = get_chunks(non_chimere[elt], chunk_size)
                perc_identity_matrix = [[]]
                for j, chunk in enumerate(segments):
                    align = nw.global_align(chunk, s[j])
                    identity =  get_identity(align)
                    perc_identity_matrix[j].append(identity)
        if not detect_chimera(perc_identity_matrix):
            non_chimere.append(squence)
            kmer_dict = get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size)
            id_seq += 1
            yield [sequence, count]

""" Regroupement glouton """

# Retourne une liste d’OTU, cette liste indiquera pour chaque séquence son occurrence (count)

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):

    otu_list = []
    sqs_count = chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)
    for elt in sqs_count:
        for seq in sqs_count:
            if (elt != seq):
                identity = get_identity([elt[0], seq[0]])
                if identity < 97:
                    otu_list.append([elt[0], elt[1]])
    return otu_list

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

# Prend une liste d’OTU et le chemin vers un fichier de sortie et affiche les OTU au format:
#>OTU_{numéro partant de 1} occurrence:{nombre d’occurrence à la déréplication}
#{séquence au format fasta}

def write_OTU(OTU_list, output_file):

    with open(output_file, "w+") as f:
        count = 1
        for i, (sequence, occurence) in enumerate(OTU_list):
            f.write(f">OTU_{count} occurrence:{occurence}\n")
            f.write(f"{fill(sequence)}")
            f.write("\n")
            count += 1

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    amplicon_file = args.amplicon_file
    if isfile(amplicon_file):
        otu_list = abundance_greedy_clustering(amplicon_file, args.minseqlen,
                  args.mincount, args.chunk_size, args.kmer_size)
        write_OTU(otu_list, args.output_file)

if __name__ == '__main__':
    main()
