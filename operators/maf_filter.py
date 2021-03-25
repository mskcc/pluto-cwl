#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run the maf_filter CWL
"""
from pluto.tools import CWLFile
from .input import generate_input
from .classes import Operator

class MafFilter(Operator):
    cwl_file = CWLFile('maf_filter.cwl')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.input = generate_input(
            self.args,
            File_keys = ['maf_file'],
            bool_keys =['is_impact', 'keep_rejects']
        )
