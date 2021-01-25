#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run TMB workflow based on CLI inputs
"""
import os
import sys
import json
from pluto.tools import TableReader, CWLFile, CWLRunner
from .input import generate_input

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

def main(
    data_clinical_file,
    assay_coverage,
    pairs_file,
    dir = None,
    verbose = True,
    func = None):
    data_clinical_file = os.path.abspath(data_clinical_file)
    pairs_file = os.path.abspath(pairs_file)

    args = {
    'data_clinical_file': data_clinical_file,
    'assay_coverage': assay_coverage,
    'pairs_file': pairs_file
    }

    input = generate_input(
        args,
        pair_template = pair_template,
        File_keys = ['data_clinical_file']
    )
    runner = CWLRunner(cwl_file = cwl_file, input = input, dir = dir, verbose = verbose)
    output_json, output_dir, output_json_file = runner.run()
    return(output_json, output_dir, output_json_file)
