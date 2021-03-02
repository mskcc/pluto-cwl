#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run the maf2vcf CWL
"""
from pluto.tools import CWLFile
from .input import generate_input
from .classes import Operator

class Maf2Vcf(Operator):
    cwl_file = CWLFile('maf2vcf.cwl')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.input = generate_input(
            self.args,
            File_keys = ['maf_file', 'ref_fasta'], File_keys_abspath = True
        )
