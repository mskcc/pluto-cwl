#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run the consensus_maf CWL
"""
from pluto.tools import CWLFile
from .input import generate_input
from .classes import Operator

class ConsensusMaf(Operator):
    cwl_file = CWLFile('consensus_maf.cwl')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        if self.args.get('maf_files_list', False):
            self.args['maf_files'] = self.args.pop('maf_files_list')
            self.input = generate_input(
                self.args,
                array_File_keys = ['maf_files']
            )
        else:
            self.input = generate_input(
                self.args,
                list_File_keys = ['maf_files']
            )
