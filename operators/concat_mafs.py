#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run the concat-mafs CWL
"""
from pluto.tools import CWLFile
from .input import generate_input
from .classes import Operator

class ConcatMafs(Operator):
    cwl_file = CWLFile('concat-mafs.cwl')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        if self.args.get('input_files_list', False):
            self.args['input_files'] = self.args.pop('input_files_list')
            self.input = generate_input(
                self.args,
                array_File_keys = ['input_files']
            )
        else:
            self.input = generate_input(
                self.args,
                list_File_keys = ['input_files']
            )
