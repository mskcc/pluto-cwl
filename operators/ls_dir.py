#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run the ls_dir cwl
"""
from pluto.tools import CWLFile
from .input import generate_input
from .classes import Operator

class LsDirCWL(Operator):
    cwl_file = CWLFile('ls_dir.cwl')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.input = {}

        if 'input_files' in self.args:
            if len(self.args['input_files']) > 0:
                self.input = generate_input(
                    self.args,
                    list_File_keys = ['input_files']
                )
