# Phylogenetic Tree Generator
#### Video Demo: https://youtu.be/ijJ4MBoS_mM 
#### Description:

This programme was created as a CS50 Final Project.

This project implements the Unweighted Pair Group Method with Arithmetic Mean (UPGMA) algorithm to build phylogenetic trees from genetic sequence data. Designed for bioinformatics applications, it provides insights into evolutionary relationships between species, strain divergence in disease outbreaks, or functionally conserved protein domains. The tool processes aligned DNA/RNA/protein sequences in FASTA format, computes pairwise genetic distances, and outputs a Newick-format tree for visualization.

Complex programmes with similar functionalities have existed for some time already. This programme is easy to use because it requires minimal input and generates simple visual data. This means it could also be used by biology educators teaching phylogenetics, students conducting introductory evolutionary analyses, or researchers needing rapid preliminary insights without deploying enterprise-scale tools.

Command line input was used for ease of use, automation and testability. FASTA files are a common format for genetic sequences as well as amino acid sequences, resulting in broad applicability of the programme; it handles aligned sequence inputs for genes (such as 16S rRNA, COI), proteins, or conserved domains. The only requirement is that the sequences be pre-aligned, which can be done using publicly available bio-informatics tools (e.g. BLAST). NumPy optimized matrix operations were applied to ensure scalability.

**Quick Start**
1. Install requirements: pip install -r requirements.txt ;
2. Run with test data: python main.py test-data/test1.fna test-data/test2.fna test-data/test3.fna ;
3. View output dendrogram: results_dendrogram.pdf ;
4. Run your own files via Command Line Interface: python main.py <input.fasta> .

**Install Dependencies**
- pip install biopython numpy scipy matplotlib

**Valid Input Files Sequences**
- Sequences must be pre-aligned in FASTA format;
- If unaligned: use alignment tools like BLAST for preparation;
- File extension must be FASTA extension

**Key Functions:**

read_seq(fasta_path):
- uses command line input of files from main();
- uses Biopython's SeqIO to open FASTA files and extract the sequence data;
- raises TypeError if a file without a FASTA extension is entered.

calculate_distances(sequences):
- returns a distance matrix comparing all original observations using pairwise-comparison output of hamming_formula(), where:
	- remaining refers to the amount of remaining iterations (remaining) == n_sequences - ref_seq;
	- ref_seq: refers to sequence in the row being iterated over a.ka. the index number of current REFERENCE sequence (sequence_1 in hamming_distance());
	- query_seq: refers to the sequence in the column of the current iteration, a.k.a. the index number of the sequence to be compared to the reference sequence;

hamming_formula(sequence1, sequence2):
- pairwise comparison of sequences;
- raises ValueError if the sequences are not of equal length (as is the case for unaligned sequences);
- outputs hamming distance between sequences.

clustering(condensed_matrix):
- uses SciPy's linkage method to create UPGMA clusters from the condensed distance matrix (contains only the upper right triangle of the matrix to avoid double comparisons).

output(cluster, labels):
- takes cluster and label input from main
- labels are the file names
- outputs a PDF of a dendrogram plot using matplotlib and the clustered data.
- the length of the two legs of a U-link in the dendrogram represents the distance between the child clusters ( [from scipy docs](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html) ) .

**Unit Testing**

Test functions were included to ensure proper error catching occurred, and the mathematical operations in hamming_formula() and calculate_distances() are accurate.
Simulated test data files test1.fna, test2.fna, test3.fna, test4.fna and test4.txt (in test-data folder) are included to test the programme functionality.

**Applications**

Although more complex, powerful tools similar to this programme are available, this programme is functional and has excellent useability for those without in-depth bioinformatics background. It can be used by students and educators to simulate or perform analysis of outbreaks like SARS-CoV-2, study conservation biology or gene family evolution, etc.. Furthermore, the scalability of the programme allows for high-level comparisons.

**Limitations**

- UPGMA complexity puts some constraint on scalability;
- Requires pre-aligned sequences;
- Does not display numerical distances in phylogenetic tree pdf. However, individual distances can be accessed by implementing a print(cluster) function in main.