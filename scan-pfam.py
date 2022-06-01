#!/usr/bin/python

import sys
import csv
import json
import requests
import collections
from os.path import exists
from tqdm import tqdm

from Bio import SeqIO

# There are many possible ways of fetching data from hmmscan, i.a.
# 1. Batch search. Send file with all sequences and asynchronously
#    wait for each to complete and download result separately.
# 2. Individual search with tsv (csv) output. Requires 2 requests,
#    one for job submit, the other for fetching tsv file.
# 3. Individual search with json output. Requires the least amount of
#    parsing and only one request to API. This is the chosen option.

def main(argv):
    if len(argv) != 2:
        print("scan_pfam.py <input_file> <output_file>")
        sys.exit()

    input_path  = argv[0]
    output_path = argv[1]

    if not exists(input_path):
        print("Input file not found.")
        sys.exit(1)

    # Scan PFAM database for domains of input sequences
    all_domains, seqs_domains = scan_pfam(input_path)

    # Write hits to csv file
    with open(output_path, 'w') as csv_file:
        filewriter = csv.writer(csv_file, delimiter=',')
        filewriter.writerow(['protein_id'] + list(all_domains))
        for seq_id, seq_domains in seqs_domains.items():
            was_seq_hit = lambda domain : 1 if domain in seq_domains else 0
            row = [seq_id] + [was_seq_hit(domain) for domain in all_domains]
            filewriter.writerow(row)

def scan_pfam(input_path):
    def fetch_hmmscan(seq_string):
        url = 'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan'
        payload = {'hmmdb': 'pfam', 'seq': seq_string}
        headers = {'Expect': '', 'Accept': 'application/json'}
        return requests.post(url, data=json.dumps(payload), headers=headers)

    all_domains = set()
    seqs_domains = collections.OrderedDict()

    print("Fetching data from hmmscan:")
    for seq_record in tqdm(list(SeqIO.parse(input_path, "fasta"))):
        response = fetch_hmmscan(str(seq_record.seq).replace("-", ""))
        hits = json.loads(response.text)['results']['hits']
        seq_domains = set()
        for hit in hits:
            for domain in hit['domains']:
                domain_id = domain['alihmmacc']
                all_domains.add(domain_id)
                seq_domains.add(domain_id)
        seqs_domains[seq_record.id] = seq_domains
    return all_domains, seqs_domains

if __name__ == "__main__":
   main(sys.argv[1:])
