#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile
sys.path.pop(0)

class TestCalcTMB(PlutoTestCase):
    cwl_file = CWLFile('calc-tmb.cwl')

    def test_calc_tmb_from_file(self):
        """
        Test case for calculating TMB by reading the number of variants from a file
        NOTE: variant entries are not parsed, only counted, so their contents do not matter here
        """
        maf_lines = [
            ['# comment 1'],
            ['# comment 2'],
            ['Hugo_Symbol', 'Chromosome'],
            ['SUFU', '1'],
            ['SUFU', '1'],
            ['SUFU', '1'],
            ['SUFU', '1'],
            ['GOT1', '2']
        ]
        input_maf_file = self.write_table(tmpdir = self.tmpdir, filename = 'input.maf', lines = maf_lines)
        self.input = {
            "input_file": {
                  "class": "File",
                  "path": input_maf_file
                },
            "output_filename":  'output.txt',
            "genome_coverage": "1000",
            "normal_id": "Sample1-N"
            }
        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'output.txt'),
                'basename': 'output.txt',
                'class': 'File',
                'checksum': 'sha1$3e5448b95b5b3e382d0ab22e0298ffd2cd20e8b9',
                'size': 7,
                'path':  os.path.join(output_dir,'output.txt')
                }
            }
        self.assertDictEqual(output_json, expected_output)

        output_file = expected_output['output_file']['path']
        with open(output_file) as fin:
            result = next(fin).strip()
        expected_result = '5000.0'
        self.assertEqual(result, expected_result)

    def test_calc_tmb_poolednormal1(self):
        """
        Test case for skipping TMB calculation if the normal is a pooled normal
        """
        maf_lines = [
            ['# comment 1'],
            ['# comment 2'],
            ['Hugo_Symbol', 'Chromosome'],
            ['SUFU', '1'],
            ['SUFU', '1'],
            ['SUFU', '1'],
            ['SUFU', '1'],
            ['GOT1', '2']
        ]
        input_maf_file = self.write_table(tmpdir = self.tmpdir, filename = 'input.maf', lines = maf_lines)
        self.input = {
            "input_file": {
                  "class": "File",
                  "path": input_maf_file
                },
            "output_filename":  'output.txt',
            "genome_coverage": "1000",
            "normal_id": "Sample1PooledNormal"
            }
        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'output.txt'),
                'basename': 'output.txt',
                'class': 'File',
                'checksum': 'sha1$7d9f637deb947080f063e9d70bdd6090968e1a7e',
                'size': 3,
                'path':  os.path.join(output_dir,'output.txt')
                }
            }
        self.assertDictEqual(output_json, expected_output)

        output_file = expected_output['output_file']['path']
        with open(output_file) as fin:
            result = next(fin).strip()
        expected_result = 'NA'
        self.assertEqual(result, expected_result)

if __name__ == "__main__":
    unittest.main()
