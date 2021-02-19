#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run the index_vcf CWL
"""
from pluto.tools import CWLFile
from .input import generate_input
from .classes import Operator

class IndexVcf(Operator):
    cwl_file = CWLFile('index_vcf.cwl')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.input = generate_input(
            self.args,
            File_keys = ['input_file'],
            File_keys_abspath = True
        )
