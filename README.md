Set of command-line tools operating on [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files suitable for:
- ```extend.py``` - protein sequence enrichment
- ```scan-pfam.py``` - protein domain finding
- ```fisher-test.py``` - measuring of protein enrichment significance


Pipfile and Pipfile.lock (for [Pipenv](https://pipenv.pypa.io/en/stable/)) with all required libraries is included. \
In addition, test script ```test.sh```, sample input data ```example.fasta``` as well as
all generated output files for this input are included.

The APIs used here are:
- [NCBI BLAST API](https://blast.ncbi.nlm.nih.gov/Blast.cgi) - for finding 
  corresponding extended protein sequences
- [HMMER hmmscan API](https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan) - for
  finding [PFAM](https://en.wikipedia.org/wiki/Pfam) protein domains

## Usage scenario
Suppose that we have a given set of protein subsequences (e.g. from 
[protein microarray](https://en.wikipedia.org/wiki/Protein_microarray))
that are known to be derived from a larger set of protein sequences.
And suppose that we are interested in what this larger set looks like, 
in particular what protein domains are present in it.

One possible solution to this problem is extending these subsequences to sequences 
that are known and very similar to them. Then a search for protein domains in these
longer sequences need to be done. In the end, we can calculate significance of this
set of domains in comparison to domains from input subsequences, using 
[Fisher's exact test](https://en.wikipedia.org/wiki/Fisher%27s_exact_test).
It will be measured as a deviation from a null hypothesis, that in both sets of
sequences, all types of domains should occur equally often.

The list of steps needed to obtain such a result is as follows:
1. You have an input FASTA file ```example.fasta``` with incomplete protein sequences.
2. Pass it to ```extend.py -i example.fasta``` to get extended protein sequences ```example-out.fasta```
from *BLAST nr Database*.
3. Pass the result file ```example-out.fasta``` to ```scan-pfam.py example-out.fasta example-out.csv```.
4. Pass the input file ```example.fasta``` to ```scan-pfam.py example.fasta example.csv```.
5. Pass both ```example.csv``` and ```example-out.csv``` to ```fisher-test.py```.
6. At the end, we get probability of this domain distribution. 

This whole set of instructions can be executed from ```test.sh example.fasta```.

## extend.py
Usage: ```extend.py -i <input_file> --evalue <evalue> --minident <minimal_percent_of_identity>```

Search in *NCBI nr (Non-redundant protein sequences)* database for protein sequences that match with at least one 
sequence from ```input_file```. 
The resulting file does not contain duplicates.
The match criterion is determined by parameters:
- minident - Minimal percent of amino acids pairs identity with respect to alignment length. Defaults to 90%.
- evalue - Maximal expected number of chance matches in a random model. Defaults to 10E-9.

## scan-pfam.py
Usage: ```scan_pfam.py <input_file> <output_file>```

Generates *.csv* file of all domains found in input file according to PFAM
database of protein domains. \
The *.csv* output file has the following structure:
- Protein sequence ids in table rows
- Protein domain ids in table columns
- At the intersection of row and column: 
    - 1 if this domain was found in this sequence
    - 0 otherwise
  
## fisher-test.py
Usage: ```fisher-test.py <file1> <file2>```

Calculate Fisher exact test for protein domain distributions from *.csv* files: ```file1``` and ```file2```.\
This is done using *fisher_exact* method from *scipy.stats* library, as well as own implementation.
