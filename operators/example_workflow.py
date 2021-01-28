#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run the example workflow
"""
import os
import copy
from pluto.tools import CWLFile
from .classes import Operator

class ExampleWorkflow(Operator):
    cwl_file = CWLFile('example_workflow.cwl')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.generate_input_data()

    def generate_input_data(self):
        input_args = copy.deepcopy(self.args)
        self.input = {
            'value': input_args['value'],
            'samples': [
            ]
        }
        for sample in input_args['sampleIDs']:
            self.input['samples'].append(
                { 'sample_id': sample }
            )
