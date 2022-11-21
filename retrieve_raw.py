import pridepy
import pandas as pd

accession = 'PXD008647'
output_folder = '/home/ir332/rds/hpc-work/hendy_ceramics_raw/'


raw_files = pridepy.Files()
resp_all_raw = raw_files.get_all_raw_file_list(accession)


def filter_blanks(resp, keep_blanks):
    new_json = []
    for file in resp:
        if ('blank' in file['fileName'] or 'wash' in file['fileName']) == keep_blanks:
            new_json.append(file)
    return new_json


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

hendy_metadata = pd.read_csv('/home/ir332/rds/hpc-work/hendy_ceramics_raw/hendy_metadata.csv')


hendy_metadata = pd.merge(
    pd.DataFrame({'File': sample_file_names}), hendy_metadata,
    on='File', how='outer', indicator=True)
hendy_metadata.sort_values('File')

for file in sample_file_names:
    print('Downloading {}...'.format(file))
    raw_files.download_file_from_ftp_by_name(
        accession=accession, file_name=file,
        output_folder=output_folder)
    print('Complete.\n')
