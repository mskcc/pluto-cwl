#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run the vcf_sort CWL
"""
from pluto.tools import CWLFile
from .input import generate_input
from .classes import Operator

class VcfSort(Operator):
    cwl_file = CWLFile('vcf_sort.cwl')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.input = generate_input(
            self.args,
            File_keys = ['vcf_file'],
            File_keys_abspath = True
        )
