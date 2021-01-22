#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the add_is_in_impact.cwl
"""
import os
import sys
import unittest
from tempfile import TemporaryDirectory

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import load_mutations, run_cwl, write_table, CWLFile
from pluto.settings import DATA_SETS, IMPACT_FILE
sys.path.pop(0)

cwl_file = CWLFile('add_is_in_impact.cwl')

class TestAddImpactCWL(unittest.TestCase):
    def test_add_impact_0(self):
        """
        Test IMPACT CWL with tiny dataset
        """
        maf_lines = [
            ['# comment 1'],
            ['# comment 2'],
            ['Hugo_Symbol'],
            ['SUFU'],
            ['GOT1'],
            ['BRCA']
        ]

        impact_lines = [
        ['BRCA', 'IMPACT468'],
        ['SUFU', 'IMPACT468'],
        ['SUFU', 'IMPACT505']
        ]

        with TemporaryDirectory() as tmpdir:
            input_maf = write_table(tmpdir = tmpdir, filename = 'input.maf', lines = maf_lines)
            impact_file = write_table(tmpdir = tmpdir, filename = 'impact.txt', lines = impact_lines)
            input_json = {
                "input_file": {
                      "class": "File",
                      "path": input_maf
                    },
                "output_filename":  'output.maf',
                "IMPACT_file": {
                      "class": "File",
                      "path": impact_file
                    },
                }
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file
                )

            expected_output = {
                'IMPACT_col_added_file': {
                    'location': 'file://' + os.path.join(output_dir, 'output.maf'),
                    'basename': 'output.maf',
                    'class': 'File',
                    'checksum': 'sha1$5c61f3977dad29ebc74966e8fc40a0278f9aab12',
                    'size': 126,
                    'path': os.path.join(output_dir, 'output.maf')
                    }
                }
            self.assertDictEqual(output_json, expected_output)

            comments, mutations = load_mutations(output_json['IMPACT_col_added_file']['path'])

            expected_comments = ['# comment 1', '# comment 2']
            self.assertEqual(comments, expected_comments)

            expected_mutations = [
                {'Hugo_Symbol': 'SUFU', 'is_in_impact': 'True', 'impact_assays': 'IMPACT468,IMPACT505'},
                {'Hugo_Symbol': 'GOT1', 'is_in_impact': 'False', 'impact_assays': '.'},
                {'Hugo_Symbol': 'BRCA', 'is_in_impact': 'True', 'impact_assays': 'IMPACT468'}
                ]
            self.assertEqual(mutations, expected_mutations)

    def test_add_impact_1(self):
        """
        Test that a maf file with is_in_IMPACT column comes out as expected
        """
        self.maxDiff = None
        input_maf = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf")

        with TemporaryDirectory() as tmpdir:
            input_json = {
                "input_file": {
                      "class": "File",
                      "path": input_maf
                    },
                "output_filename":  'output.maf',
                "IMPACT_file": {
                      "class": "File",
                      "path": IMPACT_FILE
                    },
            }

            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file
                )

            expected_output = {
                'IMPACT_col_added_file': {
                    'location': 'file://' + os.path.join(output_dir, 'output.maf'),
                    'basename': 'output.maf',
                    'class': 'File',
                    'checksum': 'sha1$1397fade2f877c2bcfca791407e328c5c48e6ff0',
                    'size': 15629589,
                    'path': os.path.join(output_dir, 'output.maf')
                    }
                }
            self.assertDictEqual(output_json, expected_output)

            # validate output mutation file contents
            with open(output_json['IMPACT_col_added_file']['path']) as fin:
                output_maf_lines = len(fin.readlines())
            self.assertEqual(output_maf_lines, 12518)

            input_comments,  input_mutations  = load_mutations(input_maf)
            output_comments, output_mutations = load_mutations(output_json['IMPACT_col_added_file']['path'])

            true_count=[row['is_in_impact'] for row in output_mutations].count('True')
            false_count=[row['is_in_impact'] for row in output_mutations].count('False')

            self.assertTrue(true_count == 8367)
            self.assertTrue(false_count == 4147)

            # check that its got two extra columns in the output
            self.assertTrue(len(input_mutations[1])+2==len(output_mutations[1]))

if __name__ == "__main__":
    unittest.main()
