from Bio import SeqIO
import sys
import numpy as np


def main():
    files = sys.argv[1:]
    sequence = []
    for file in files:
        sequence += [str(read_seq(file))]
    calculate_distance(sequence)


def read_seq(fasta_path):
    # opens a FASTA file and returns sequence
    with open(fasta_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            return record.seq


def hamming_formula(sequence_1, sequence_2):  # pair-wise comparison for hamming distance
    if len(sequence_1) != len(sequence_2):
        raise ValueError("Sequences must be of equal length.")
    d_min = 0
    for n in range(len(sequence_1)):
        if sequence_1[n] != sequence_2[n]:
            d_min += 1
    return d_min


def calculate_distance(sequences):
    n_sequences = len(sequences)
    distances_matrix = np.zeros((n_sequences, n_sequences))  # creates numpy array filled with zeros of length/width = n_sequences

    for ref_seq in range(n_sequences - 1):  # iterates over sequences until n-1

        # print("Sequence Index: " + str(ref_seq))  # indicates index of current REFERENCE sequence (sequence_1 in hamming_distance())
        for remaining in range(1, n_sequences - ref_seq):  # the amount of remaining iterations (remaining) == n_sequences - ref_seq
            query_seq = ref_seq + remaining  # query_seq is the index number of the sequence to be compared to the reference sequence
            # print("Query Sequence: " + str(query_seq))  # indicates index of current query sequence (sequence_2 in hamming_distance())
            hamming_distance = hamming_formula(sequences[ref_seq], sequences[query_seq])

            distances_matrix[ref_seq][query_seq] = (hamming_distance / len(sequences[ref_seq]))  # upper right triangle #normalization via /len(sequences[ref_seq]
            distances_matrix[query_seq][ref_seq] = (hamming_distance / len(sequences[ref_seq]))  # bottom left triangle (inverse)

    matrix_relations(distances_matrix)  # calculate_distances() returns distances_matrix, inserts distances_matrix value into matrix_relations parameter


def matrix_relations(matrix):
    size = len(matrix)
    matrix_indices = np.triu_indices(size, 1) # return 2 arrays with x and y coordinates
    for index in range(len(matrix_indices[0])): # run this loop as many times as there are items in the x list
        x = matrix_indices[0][index] # [0] --> the x array ; [index] --> specific x value
        y = matrix_indices[1][index] # [1] --> the y array ; [index] --> specific y value
        value_matrix = matrix[x][y]
        print(f'x:{x} Y:{y}')
        print(matrix[x][y])
    # for x in matrix_indices[0]:
    #     for y in matrix_indices[1]:
    #         print(f'x:{x} Y:{y}')
    #         print(matrix[x][y])

    #triu_matrix = np.triu(matrix) #returns array with only upper right triangle
    # https: // numpy.org / doc / stable / reference / generated / numpy.triu.html  # numpy.triu
    #matrix = np.triu_indices
    # https: // numpy.org / doc / stable / reference / generated / numpy.triu_indices.html --> Upper Triangle
    #min_value = np.min(triu_matrix)
    # print(matrix_indices)
    # print(matrix)
    # print(matrix_changed)



if __name__ == "__main__":
    main()
