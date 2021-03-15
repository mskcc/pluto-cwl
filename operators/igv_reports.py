#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run the igv-reports CWL
"""
from pluto.tools import CWLFile
from .input import generate_input
from .classes import Operator

class IgvReports(Operator):
    cwl_file = CWLFile('igv-reports.cwl')

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        if self.args['bam_files'] is None:
            self.args['bam_files'] = []
        if self.args['vcf_gz_files'] is None:
            self.args['vcf_gz_files'] = []

        self.input = generate_input(
            self.args,
            list_File_keys = ['vcf_gz_files', 'bam_files'],
            File_keys = ['ref_fasta', 'sites'],
            File_keys_abspath = True
        )
