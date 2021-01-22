#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for running the CWL workflows
"""
import csv
import argparse
from pluto.tools import write_table, dicts2lines
from operators import tmb_workflow

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

def main():
    """
    Main function for CLI parsing

    $ ./run.py tmb_workflow --data-clinical examples/data_clinical.txt --assay-coverage 10000 --pairs examples/pairs.tsv
    """
    parser = argparse.ArgumentParser(description = '')

    subparsers = parser.add_subparsers(help ='Sub-commands available', required = True)

    _tmb_workflow = subparsers.add_parser('tmb_workflow', help = 'Help goes here')
    _tmb_workflow.add_argument('--data-clinical', dest = 'data_clinical_file', required = True, help = '')
    _tmb_workflow.add_argument('--assay-coverage', dest = 'assay_coverage', required = True, help = '')
    _tmb_workflow.add_argument('--pairs', dest = 'pairs_file', required = True, help = '')
    _tmb_workflow.set_defaults(func = tmb_workflow.main)

    _generate_data_clinical = subparsers.add_parser('data_clinical_template', help = 'Help goes here')
    _generate_data_clinical.add_argument('--output', dest = 'output_file', default = "data_clinical.txt", help = '')
    _generate_data_clinical.set_defaults(func = generate_data_clinical)

    _generate_pairs = subparsers.add_parser('pairs_template', help = 'Help goes here')
    _generate_pairs.add_argument('--output', dest = 'output_file', default = "pairs.tsv", help = '')
    _generate_pairs.set_defaults(func = generate_pairs)

    args = parser.parse_args()
    args.func(**vars(args))

if __name__ == '__main__':
    main()
