#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run the snp-pileup-wrapper CWL
"""
from pluto.tools import CWLFile
from .input import generate_input
from .classes import Operator

class SnpPileup(Operator):
    cwl_file = CWLFile('snp-pileup-wrapper.cwl')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.input = generate_input(
            self.args,
            File_keys = ['snps_vcf', 'normal_bam', 'tumor_bam'],
            File_keys_abspath = True
        )
