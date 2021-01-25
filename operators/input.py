"""
Module for converting Python dict into CWL input JSON dict
"""
from pluto.tools import TableReader, write_table, dicts2lines
import copy

def generate_pairs(pairs_file, pair_template):
    """
    Parse a pairs samplesheet file to create a list of pair records for CWL 'pairs' input
    """
    table_reader = TableReader(pairs_file)
    comments = table_reader.comment_lines
    fieldnames = table_reader.get_fieldnames()
    records = [ rec for rec in table_reader.read() ]

    pairs = []
    for record in records:
        # start wit a copy of the the template
        pair = copy.deepcopy(pair_template)
        # add the pair_maf entry if it was included
        if 'pair_maf' in pair_template:
            pair['pair_maf']['path'] = record.pop('pair_maf')
        # add the snp_pileup if it was included
        if 'snp_pileup' in pair_template:
            pair['snp_pileup']['path'] = record.pop('snp_pileup')
        # add all the other items from the input record
        for k, v in record.items():
            if k in pair:
                pair[k] = v
        pairs.append(pair)
    return(pairs)

def generate_input(
    args,
    bool_keys = None,
    File_keys = None,
    array_File_keys = None,
    pair_template = None # requires 'pairs_file' in args
    ):
    """
    Build the input JSON from the args passed
    All required keys should be present in the input args
    """
    if bool_keys is None:
        bool_keys = []
    if File_keys is None:
        File_keys = []
    if array_File_keys is None:
        array_File_keys = []
    if bool_keys is None:
        bool_keys = []

    input = copy.deepcopy(args)

    # convert input args into the correct output formats
    for key in bool_keys:
        input[key] = bool(input[key])

    for key in File_keys:
        path = input[key]
        input[key] = {'class':'File', 'path': path}

    # input value is a single txt file that contains one filepath per line to build into array of File
    for key in array_File_keys:
        path = input[key]
        input[key] = []
        lines = []
        with open(path) as f:
            lines = [ line.strip() for line in f ]
        for line in lines:
            d = {'class':'File', 'path': line}
            input[key].append(d)

    # fill in the 'pairs' field with the values parsed from the pairs_file samplesheet
    if pair_template:
        pairs_file = input.pop('pairs_file')
        input['pairs'] = generate_pairs(pairs_file, pair_template)

    return(input)


def generate_data_clinical(func = None, output_file = "data_clinical.txt"):
    """
    create a blank base data clinical file template
    """
    data_clinical_lines = [
    ['#SAMPLE_ID'],
    ['#SAMPLE_ID'],
    ['#STRING'],
    ['#1'],
    ['SAMPLE_ID']
    ]
    write_table('', '', lines = data_clinical_lines, filepath = output_file)

def generate_pairs(func = None, output_file = "pairs.tsv"):
    """
    create a blank template pairs file for use with pipelines
    """
    lines = [[
    'tumor_id',
    'normal_id',
    'pair_id',
    'pair_maf',
    'snp_pileup'
    ]]
    write_table('', '', lines = lines, filepath = output_file)

def generate_sample_summary(func = None, output_file = "sample_summary.txt"):
    lines = [[
        'Auto-status',
        'Sample',
        'Unexpected Match(es)',
        'Unexpected Mismatch(es)',
        'Major Contamination',
        'Minor Contamination',
        'Coverage',
        'Duplication',
        'Library Size (millions)',
        'On Bait Bases (millions)',
        'Aligned Reads (millions)',
        'Insert Size Peak'
    ]]
    write_table('', '', lines = lines, filepath = output_file)
