#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the add_is_in_impact.cwl
"""
import os
import sys
import unittest
from datasets import (
    DATA_SETS
)
from pluto import (
    PlutoTestCase, 
    CWLFile,
    IMPACT_FILE,
    OFile
)



class TestAddIsInImpact(PlutoTestCase):
    cwl_file = CWLFile('add_is_in_impact.cwl')

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

        input_maf = self.write_table(tmpdir = self.tmpdir, filename = 'input.maf', lines = maf_lines)
        impact_file = self.write_table(tmpdir = self.tmpdir, filename = 'impact.txt', lines = impact_lines)
        self.input = {
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

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'IMPACT_col_added_file': OFile(name = 'output.maf', hash = '5c61f3977dad29ebc74966e8fc40a0278f9aab12', size = 126, dir = output_dir)
            }
        self.assertCWLDictEqual(output_json, expected_output)

        self.assertMutFileContains(
            filepath = expected_output['IMPACT_col_added_file']['path'],
            expected_comments = ['# comment 1', '# comment 2'],
            expected_mutations = [
                {'Hugo_Symbol': 'SUFU', 'is_in_impact': 'True', 'impact_assays': 'IMPACT468,IMPACT505'},
                {'Hugo_Symbol': 'GOT1', 'is_in_impact': 'False', 'impact_assays': '.'},
                {'Hugo_Symbol': 'BRCA', 'is_in_impact': 'True', 'impact_assays': 'IMPACT468'}
                ],
            identical = True
        )

    def test_add_impact_1(self):
        """
        Test that a maf file with is_in_IMPACT column comes out as expected
        """
        self.maxDiff = None
        input_maf = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf")

        self.input = {
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

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'IMPACT_col_added_file': OFile(name = "output.maf", size = 15629589, hash = "1397fade2f877c2bcfca791407e328c5c48e6ff0", dir = output_dir)
            }
        self.assertCWLDictEqual(output_json, expected_output)

        # validate output mutation file contents
        with open(expected_output['IMPACT_col_added_file']['path']) as fin:
            output_maf_lines = len(fin.readlines())
        self.assertEqual(output_maf_lines, 12518)

        input_comments,  input_mutations  = self.load_mutations(input_maf)
        output_comments, output_mutations = self.load_mutations(expected_output['IMPACT_col_added_file']['path'])

        true_count=[row['is_in_impact'] for row in output_mutations].count('True')
        false_count=[row['is_in_impact'] for row in output_mutations].count('False')

        self.assertTrue(true_count == 8367)
        self.assertTrue(false_count == 4147)

        # check that its got two extra columns in the output
        self.assertTrue(len(input_mutations[1]) + 2 ==len(output_mutations[1]))

