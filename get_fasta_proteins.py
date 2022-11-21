import requests
import sys
import os
import re
import argparse
from Bio import SeqIO
import json
from io import StringIO
from itertools import chain


def divide_chunks(li, n):
    # looping till length l
    for i in range(0, len(li), n):
        yield li[i:i + n]


def map_orgnames(species):

    tax_params = {
        'pageNumber': 1,
        'pageSize': 1,
        'searchType': 'EQUALSTO',
        'fieldName': 'SCIENTIFICNAME',

    }

    taxid_list = []
    for sp in species:
        requestURL = ("https://www.ebi.ac.uk/proteins/api/taxonomy/name/" +
                      "{}/node?".format(sp))
        r = requests.get(requestURL, params=tax_params,
                         headers={"Accept": "application/json"})
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        result_json = json.loads(r.text)

        taxid = result_json['taxonomies'][0]['taxonomyId']
        print('{}: {}'.format(sp, taxid))
        taxid_list.append(str(taxid))
    return taxid_list


def get_fasta(requestURL, pars):
    fasta_request = requests.get(requestURL, params=pars, headers={"Accept": "text/x-fasta"})
    if not fasta_request.ok:
        fasta_request.raise_for_status()
        sys.exit()

    fasta_string = fasta_request.text
    fasta_io = StringIO(fasta_string)

    seq_entries = list(SeqIO.parse(fasta_io, format='fasta'))

    return seq_entries


def get_prots_by_sp_and_name(taxid_list, genes, maxsp, maxgenes):

    taxid_chunks = divide_chunks(taxid_list, maxsp)
    genes_chunks = divide_chunks(genes, maxgenes)
    requestURL = "https://www.ebi.ac.uk/proteins/api/proteins?"
    seq_entries = []
    for tc in taxid_chunks:
        for gc in genes_chunks:
            pars = {'offset': '0',
                    'size': '100',
                    'isoform': '2',
                    'exact_gene': ','.join(gc),
                    'taxid': ','.join(tc)}
            seq_chunk = get_fasta(requestURL, pars)
            seq_entries = chain(seq_entries, seq_chunk)

    return list(seq_entries)


def get_prots_by_acc(accessions, maxacc):

    requestURL = "https://www.ebi.ac.uk/proteins/api/proteins?"
    seq_entries = []
    acc_chunks = divide_chunks(accessions, maxacc)
    for ac in acc_chunks:
        pars = {'offset': '0',
                'size': '100',
                'accession': ','.join(ac)}
        seq_chunk = get_fasta(requestURL, pars)
        seq_entries = chain(seq_entries, seq_chunk)

    return list(seq_entries)


def filter_seqs(seq_entries, regex):
    filt_seqs = []
    regex = re.compile(regex, re.IGNORECASE)
    for record in seq_entries:
        if not regex.search(record.description):
            filt_seqs.append(record)
    return filt_seqs


def main():
    arg_parser = argparse.ArgumentParser(description='Script to query Uniprot and retrieve sequences')
    cwd = os.getcwd()
    arg_parser.add_argument(
        '-s', '--species', required=False, action='store', default=None,
        help='File listing species in separate lines')
    arg_parser.add_argument(
        '-ms', '--maxsp', required=False, action='store', default=20, type=int,
        help='Maximum species in API request')
    arg_parser.add_argument(
        '-g', '--genes', required=False, action='store', default=None,
        help='File listing genes names in separate lines')
    arg_parser.add_argument(
        '-mg', '--maxgenes', required=False, action='store', default=20, type=int,
        help='Maxsimum genes in API request')
    arg_parser.add_argument(
        '-a', '--accessions', required=False, action='store', default=None,
        help='File listing protein accesions in separate lines')
    arg_parser.add_argument(
        '-ma', '--maxacc', required=False, action='store', default=100, type=int,
        help='Maximum accessions in API request')
    arg_parser.add_argument(
        '-o', '--outfile', required=False, default=cwd + '/uniprot_seqs.fasta',
        help='Outut file')
    args = arg_parser.parse_args()

    species_f = args.species
    maxsp = args.maxsp
    genes_f = args.genes
    maxgenes = args.maxgenes
    acc_f = args.accessions
    maxacc = args.maxacc
    outfile = args.outfile

    if species_f is None and genes_f is None and acc_f is None:
        sys.exit(1)

    seq_sp = []
    if species_f is not None and genes_f is not None:
        with open(species_f) as spf:
            species = spf.readlines()
        species = [s.rstrip() for s in species]

        with open(genes_f) as gf:
            genes = gf.readlines()
        genes = [g.rstrip() for g in genes]

        taxid_list = map_orgnames(species)
        seq_sp = get_prots_by_sp_and_name(taxid_list, genes, maxsp, maxgenes)
        seq_sp = filter_seqs(seq_sp, '(fragment)|(isoform)')

    seq_acc = []
    if acc_f is not None:
        with open(acc_f) as accf:
            accessions = accf.readlines()
        accessions = [a.rstrip() for a in accessions]
        seq_acc = get_prots_by_acc(accessions, maxacc)

    # Combine fastas
    seq_entries = seq_sp + seq_acc

    print('Writing FASTA file in {} ...'.format(outfile))
    SeqIO.write(seq_entries, outfile, format='fasta')


if __name__ == '__main__':
    main()
