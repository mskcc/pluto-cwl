#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run the concat-tables_dir CWL
"""
import copy
from pluto.tools import CWLFile
from .input import generate_input
from .classes import Operator

class ConcatTablesDirCWL(Operator):
    cwl_file = CWLFile('concat-tables_dir.cwl')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        if self.args.get('input_files_list', False):
            self.args['input_files'] = self.args.pop('input_files_list')
            self.input = generate_input(
                self.args,
                array_File_keys = ['input_files'],
                bool_keys = ['comments']
            )
        else:
            self.input = generate_input(
                self.args,
                list_File_keys = ['input_files'],
                bool_keys = ['comments']
            )
