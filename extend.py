#!/usr/bin/python

import pickle
import sys, getopt
from tqdm import tqdm
from os.path import exists, splitext

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML


Entrez.email = "your@email.com"

def search_ncbi(input_path, expect, min_perc_ident):
    def query_ncbi(seq_string):
        return NCBIWWW.qblast("blastp", "nr", seq_string, expect=expect,
                              perc_ident=min_perc_ident)

    filename, _ = splitext(input_path)
    output_path = filename + "-out.fasta"
    blastp_result_pkl = filename + ".pkl"

    # Get sequence alignments from NCBI API
    if not exists(blastp_result_pkl):
        blastp_result = []
        print("Fetching data from NCBI: ")
        for seq in tqdm(list(SeqIO.parse(input_path, "fasta"))):
            seq_string = seq.format("fasta")
            result_handle = query_ncbi(seq_string)
            blast_records = NCBIXML.parse(result_handle)
            blastp_result.append(list(blast_records)[0])
            result_handle.close()

            # Save current result to pickle file
            pickle_file = open(blastp_result_pkl, 'wb')
            pickle.dump(blastp_result, pickle_file)
            pickle_file.close()

    # Load sequence alignments
    file = open(blastp_result_pkl, 'rb')
    blastp_result = pickle.load(file)
    file.close()

    # Get hit sequences from alignments
    result_seqrecords = []
    result_seqs = set()

    for n, blast_record in enumerate(blastp_result):
        for alignment in blast_record.alignments:
            # Take only first hsp, because we are not interested
            # in the location, but in the entire sequence
            hsp = alignment.hsps[0]

            # Check if the matched sequence is new
            if hsp.sbjct in result_seqs:
                continue

            perc_ident = 100 * hsp.identities / hsp.align_length
            # Need to filter identity percent because
            # NCBI perc_ident parameter is deprecated
            if (perc_ident < min_perc_ident):
                continue

            hsp_data = {
                "matched": blast_record.query,
                "score": hsp.score,
                "expect": hsp.expect,
                "perc_ident": str(round(perc_ident, 2)) + "%",
            }

            def dict_to_str(d):
                return str(d).replace("{","").replace("}","").replace("'","")

            seq_rec = SeqRecord(
                Seq(hsp.sbjct),
                id = alignment.hit_id,
                description = dict_to_str(hsp_data),
            )

            result_seqs.add(hsp.sbjct)
            result_seqrecords.append(seq_rec)

    # Write hit sequences to output file
    SeqIO.write(result_seqrecords, output_path, "fasta")

def main(argv):
    help_string = "extend.py -i <input_file> --evalue <evalue>" \
                  " --minident <minimal_percent_of_identity>"
    input_path  = ""
    expect      = 10e-10 # Expected number of chance matches in a random model
    min_perc_ident = 90  # Minimum hit identity percent

    try:
        opts, args = getopt.getopt(argv, "hi:", ["evalue=", "minident="])
    except getopt.GetoptError:
        print(help_string)
        sys.exit(1)
    for opt, arg in opts:
        if opt == "-h":
            print(help_string)
            sys.exit()
        elif opt == "-i":
            input_path = arg
        elif opt == "--evalue":
            expect = arg
        elif opt == "--minident":
            min_perc_ident = arg

    if not exists(input_path):
        print(help_string)
        sys.exit(1)

    search_ncbi(input_path, expect, min_perc_ident)

if __name__ == "__main__":
   main(sys.argv[1:])
