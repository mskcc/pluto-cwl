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
        self.args['bam_files'] = []
        self.args['maf_files'] = []
        self.args['sample_ids'] = []
        samplesheet = os.path.abspath(samplesheet)
        # load the records from the samplesheet
        with open(samplesheet) as f:
            reader = csv.DictReader(f, delimiter = '\t')
            for row in reader:
                self.args['bam_files'].append(row['bam_file'])
                self.args['maf_files'].append(row['maf_file'])
                self.args['sample_ids'].append(row['sample_id'])

        # convert all inputs to CWL input format
        self.input = generate_input(
            self.args,
            File_keys = ['ref_fasta'],
            list_File_keys = ['bam_files', 'maf_files']
        )
