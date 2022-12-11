#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the maf2vcf_gz cwl
"""
import os
import sys
import unittest
from collections import OrderedDict

PARENT_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, PARENT_DIR)
from pluto import (
    PlutoTestCase, 
    CWLFile, 
    TableReader
)
sys.path.pop(0)

class TestMaf2VcfGz(PlutoTestCase):
    cwl_file = CWLFile('maf2vcf_gz_workflow.cwl')

    def test_convert_maf_to_vcf(self):
        """
        Convert a single input maf file into a .vcf.gz with a .tbi
        """
        self.maxDiff = None
        self.input = {
            "maf_file": {"class": "File", "path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample4.Sample3.muts.maf") },
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']},
        }

        output_json, output_dir = self.run_cwl()

        # do not include size and checksum since they are not consistent with .gz
        output_json['output_file'].pop('checksum')
        output_json['output_file'].pop('size')
        output_json['output_file']['secondaryFiles'][0].pop('checksum')
        output_json['output_file']['secondaryFiles'][0].pop('size')

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'variants.vcf.gz'),
                'basename': 'variants.vcf.gz',
                'class': 'File',
                # 'checksum': 'sha1$dcda65da7665683dcf97ed9c3989fd75c7a839c8',
                # 'size': 14885,
                'secondaryFiles': [
                    {
                        'basename': 'variants.vcf.gz.tbi',
                        'location': 'file://' + os.path.join(output_dir,'variants.vcf.gz.tbi'),
                        'class': 'File',
                        # 'checksum': 'sha1$c537928ca8bd7f33d19302d42c84ed6370687fca',
                        # 'size': 9320,
                        'path': os.path.join(output_dir,'variants.vcf.gz.tbi')
                    }
                ],
                'path': os.path.join(output_dir,'variants.vcf.gz') }}

        self.assertCWLDictEqual(output_json, expected_output)




