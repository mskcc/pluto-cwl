#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run the deduplicate-maf CWL
"""
from pluto.tools import CWLFile
from .input import generate_input
from .classes import Operator

class DeduplicateMaf(Operator):
    cwl_file = CWLFile('deduplicate-maf.cwl')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.input = generate_input(
            self.args,
            File_keys = ['input_file']
        )
