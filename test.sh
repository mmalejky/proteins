#!/bin/sh

input_basename="${1%.*}"

python extend.py -i "$1" &&
python scan-pfam.py "$1" "$input_basename.csv" &&
python scan-pfam.py "$input_basename-out.fasta" "$input_basename-out.csv" &&
python fisher-test.py "$input_basename.csv" "$input_basename-out.csv"
