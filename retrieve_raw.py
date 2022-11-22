import pridepy
import pandas as pd
import argparse




def filter_blanks(resp, keep_blanks):
    new_json = []
    for file in resp:
        if ('blank' in file['fileName'] or 'wash' in file['fileName']) == keep_blanks:
            new_json.append(file)
    return new_json







def main():
    arg_parser = argparse.ArgumentParser(description='Script to query Uniprot and retrieve sequences')

    arg_parser.add_argument(
        '-a', '--accession', required=True, action='store', default=None,
        help='PRIDE accession')
    arg_parser.add_argument(
        '-o', '--outfolder', required=True, action='store', default=None,
        help='Output folder where RAW files are saved')

    args = arg_parser.parse_args()

    accession = args.accession
    output_folder = args.outfolder
    raw_files = pridepy.Files()
    resp_all_raw = raw_files.get_all_raw_file_list(accession)

    raw_samples = filter_blanks(resp_all_raw, keep_blanks=False)
    raw_blank_wash = filter_blanks(resp_all_raw, keep_blanks=True)

    sample_file_names = []
    for file in raw_samples:
        sample_file_names.append(file['fileName'])

    blank_file_names = []
    for file in raw_blank_wash:
        blank_file_names.append(file['fileName'])
    blank_file_names = pd.DataFrame({'File': blank_file_names})
    blank_file_names.sort_values('File')

    for file in sample_file_names:
        print('Downloading {}...'.format(file))
        raw_files.download_file_from_ftp_by_name(
            accession=accession, file_name=file,
            output_folder=output_folder)
        print('Complete.\n')


if __name__ == '__main__':
    main()
