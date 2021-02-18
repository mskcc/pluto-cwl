#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the add_af.cwl
"""
import os
import sys
import json
import unittest
from tempfile import TemporaryDirectory

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import load_mutations, run_cwl, write_table, CWLFile

sys.path.pop(0)

cwl_file = CWLFile('add_af.cwl')

class TestAddImpactCWL(unittest.TestCase):
    def test_add_af(self):
        """
        Test IMPACT CWL with tiny dataset
        """
        maf_lines = [
            ['# comment 1'],
            ['# comment 2'],
            ['Hugo_Symbol', 't_depth', 't_alt_count'],
            ['SUFU', '100', '75'],
            ['GOT1', '100', '1'],
            ['SOX9', '100', '0'],
        ]

        with TemporaryDirectory() as tmpdir:
            input_maf = write_table(tmpdir = tmpdir, filename = 'input.maf', lines = maf_lines)
            input_json = {
                "input_file": {
                      "class": "File",
                      "path": input_maf
                    },
                "output_filename":  'output.maf',
                }
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file
                )

            expected_output = {
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'output.maf'),
                    'basename': 'output.maf',
                    'class': 'File',
                    'checksum': 'sha1$39de59ad5d736db692504012ce86d3395685112e',
                    'size': 109,
                    'path': os.path.join(output_dir, 'output.maf')
                    }
                }
            self.assertDictEqual(output_json, expected_output)

            comments, mutations = load_mutations(output_json['output_file']['path'])

            expected_comments = ['# comment 1', '# comment 2']
            self.assertEqual(comments, expected_comments)

            expected_mutations = [
                {'Hugo_Symbol': 'SUFU', 't_depth': '100', 't_alt_count':'75', 't_af': '0.75'},
                {'Hugo_Symbol': 'GOT1', 't_depth': '100', 't_alt_count':'1', 't_af': '0.01'},
                {'Hugo_Symbol': 'SOX9', 't_depth': '100', 't_alt_count':'0', 't_af': '0.0'}
                ]
            self.assertEqual(mutations, expected_mutations)

if __name__ == "__main__":
    unittest.main()
