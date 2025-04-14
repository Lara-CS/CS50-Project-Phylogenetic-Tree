from Bio import SeqIO
import sys
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from matplotlib import pyplot as plt
import os.path


def main():
    files = sys.argv[1:]
    sequence = []
    labels = []
    for file in files:
        try:
            sequence += [str(read_seq(file))]
        except TypeError:
            sys.exit("A non-FASTA file was entered for analysis.")
        except FileNotFoundError:
            sys.exit("An invalid file name was entered, or an unavailable file was entered for analysis")
        labels += [os.path.basename(file).split('.')[0]]
    distances_matrix = calculate_distance(sequence)
    condensed_matrix = condensation(distances_matrix)
    cluster = clustering(condensed_matrix)
    output(cluster, labels)


def read_seq(fasta_path):
    # opens a FASTA file and returns sequence
    if not fasta_path.endswith((".fasta", ".fas", ".fa", ".fna", ".ffn", ".faa", ".mpfa", ".frn")):
        raise TypeError()
    with open(fasta_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return record.seq


def condensation(distances_matrix):
    condensed_matrix = squareform(distances_matrix)
    return condensed_matrix


def clustering(condensed_matrix):
    cluster = linkage(condensed_matrix, method= 'average') # average --> UPGMA clustering
    return cluster


def output(cluster, labels):
    fig = plt.figure(figsize=(25, 10))
    plt.title("Phylogenetic Tree of Sequences")
    dn = dendrogram(cluster, orientation='left', labels=labels) # to be added: parameter --> labels=sequence_name
    plt.savefig('results_dendrogram.pdf')

def calculate_distance(sequences):
    n_sequences = len(sequences)
    distances_matrix = np.zeros((n_sequences, n_sequences))  # creates numpy array filled with zeros of length/width = n_sequences

    for ref_seq in range(n_sequences - 1):  # iterates over sequences until n-1

        for remaining in range(1, n_sequences - ref_seq):  # the amount of remaining iterations (remaining) == n_sequences - ref_seq
            query_seq = ref_seq + remaining  # query_seq is the index number of the sequence to be compared to the reference sequence
            hamming_distance = hamming_formula(sequences[ref_seq], sequences[query_seq])

            distances_matrix[ref_seq][query_seq] = (hamming_distance / len(sequences[ref_seq]))
            distances_matrix[query_seq][ref_seq] = (hamming_distance / len(sequences[ref_seq]))

    return distances_matrix


def hamming_formula(sequence_1, sequence_2):  # pair-wise comparison for hamming distance
    if len(sequence_1) != len(sequence_2):
        raise ValueError("Sequences must be of equal length.")
    d_min = 0
    for n in range(len(sequence_1)):
        if sequence_1[n] != sequence_2[n]:
            d_min += 1
    return d_min


if __name__ == "__main__":
    main()