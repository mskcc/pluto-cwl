#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run the samples_fillout_workflow CWL
"""
import os
import csv
from pluto.tools import CWLFile
from .input import generate_input
from .classes import Operator

class SamplesFillout(Operator):
    cwl_file = CWLFile('samples_fillout_workflow.cwl')
    def __init__(self, samplesheet, **kwargs):
        super().__init__(**kwargs)
        samplesheet = os.path.abspath(samplesheet)
        with open(samplesheet) as f:
            reader = csv.DictReader(f, delimiter = '\t')
            rows = [ row for row in reader ]

        self.input = generate_input(
            self.args,
            File_keys = ['ref_fasta']
        )
        self.input['samples'] = []
        for row in rows:
            sample = {
                'sample_id': row['sample_id'],
                'maf_file': {'class': 'File', 'path': row['maf_file']},
                'bam_file': {'class': 'File', 'path': row['bam_file']}
            }
            self.input['samples'].append(sample)
