#!/usr/bin/python

import sys
import csv
import numpy as np
from scipy.stats import fisher_exact
from scipy.special import comb

def main(argv):
    if len(argv) != 2:
        print("fisher-test.py <file1> <file2>")
        sys.exit(1)

    files_count = len(argv)
    proteins = np.zeros(files_count, dtype=int)
    hits     = np.zeros(files_count, dtype=int)

    for file_number, file_path in enumerate(argv):
        file = open(file_path, mode='r')
        csv_reader = csv.reader(file)
        next(csv_reader) # Skip first line
        for row in csv_reader:
            proteins[file_number] += 1
            hits[file_number]     += 1 if "1" in row else 0
        file.close()


    result1 = my_fisher(proteins, hits)
    print("My fisher test result:", result1)

    no_hits = proteins - hits
    contingency_table = [hits, no_hits]
    result2 = fisher_exact(contingency_table)[1]
    print("Fisher exact test result:", result2)

def my_fisher(proteins, hits):
    N = np.sum(proteins) # Total proteins
    n = np.sum(hits)     # Total hits

    # Larger sample file index
    larger_idx = 0 if proteins[0] > proteins[1] else 1
    K = proteins[larger_idx] # Proteins in larger sample file
    x = hits[larger_idx]     # Hits in larger sample file
    print("N:", N, "K:", K, "n:", n, "x:", x)

    print("file_no", "seqs", "hits", sep="\t")
    for i in range(len(proteins)):
        print(i, proteins[i], hits[i], sep="\t")

    i = np.arange(x, min(K, n)+1)
    return np.sum(comb(K, i) * comb(N-K, n-i)) / comb(N, n)

if __name__ == "__main__":
   main(sys.argv[1:])
