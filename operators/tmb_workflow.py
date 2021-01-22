#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run TMB workflow based on CLI inputs
"""
import os
import sys
import json
from pluto.tools import TableReader, CWLFile, CWLRunner

cwl_file = CWLFile('tmb_workflow.cwl')

pair_template = {
        "pair_maf": {
            "path": None,
            "class": "File"
        },
        "pair_id": None,
        "tumor_id": None,
        "normal_id": None
    }

input_template = {
    'data_clinical_file': {
          "class": "File",
          "path": None
    },
    "assay_coverage":  None,
    "pairs": [ pair_template ]
}

def generate_pairs(pairs_file):
    table_reader = TableReader(pairs_file)
    comments = table_reader.comment_lines
    fieldnames = table_reader.get_fieldnames()
    records = [ rec for rec in table_reader.read() ]

    pairs = []
    for record in records:
        pair = {**pair_template}
        pair['pair_maf']['path'] = record.pop('pair_maf')
        for k, v in record.items():
            if k in pair:
                pair[k] = v
        pairs.append(pair)
    return(pairs)

def generate_input(data_clinical_file, assay_coverage, pairs):
    input = {**input_template}
    input['pairs'] = pairs
    input['data_clinical_file']['path'] = data_clinical_file
    input['assay_coverage'] = assay_coverage
    return(input)

def main(data_clinical_file, assay_coverage, pairs_file, func = None):
    data_clinical_file = os.path.abspath(data_clinical_file)
    pairs_file = os.path.abspath(pairs_file)
    pairs = generate_pairs(pairs_file)
    input = generate_input(data_clinical_file, assay_coverage, pairs)
    runner = CWLRunner(cwl_file, input)
    runner.run()
