#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run TMB workflow based on CLI inputs
"""
import os
from pluto.tools import CWLFile
from .input import generate_input
from .classes import Operator

class TMBWorkflow(Operator):
    pair_template = {
            "pair_maf": {
                "path": None,
                "class": "File"
            },
            "pair_id": None,
            "tumor_id": None,
            "normal_id": None
        }
    cwl_file = CWLFile('tmb_workflow.cwl')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.parse_kwargs(**kwargs)
        self.generate_input_data()

    def parse_kwargs(self, **kwargs):
        super().parse_kwargs(**kwargs)
        self.data_clinical_file = os.path.abspath(kwargs.pop('data_clinical_file'))
        self.assay_coverage = kwargs.pop('assay_coverage')
        self.pairs_file = os.path.abspath(kwargs.pop('pairs_file'))

    def generate_input_data(self):
        args = {
        'data_clinical_file': self.data_clinical_file,
        'assay_coverage': self.assay_coverage,
        'pairs_file': self.pairs_file
        }
        self.input = generate_input(
            args,
            pair_template = self.pair_template,
            File_keys = ['data_clinical_file']
        )
