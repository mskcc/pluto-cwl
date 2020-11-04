#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the maf_filter.cwl
"""
import os
import json
import unittest
from tempfile import TemporaryDirectory, NamedTemporaryFile

# relative imports, from CLI and from parent project
if __name__ != "__main__":
    from .tools import run_command, load_mutations, run_cwl, write_table
    from .settings import CWL_DIR, CWL_ARGS, DATA_SETS, ARGOS_VERSION_STRING, IS_IMPACT, IMPACT_FILE

if __name__ == "__main__":
    from tools import run_command, load_mutations, run_cwl, write_table
    from settings import CWL_DIR, CWL_ARGS, DATA_SETS, ARGOS_VERSION_STRING, IS_IMPACT, IMPACT_FILE

cwl_file = os.path.join(CWL_DIR, 'add_is_in_impact.cwl')

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
            ['GOT1']
        ]

        impact_lines = [
        ['SUFU']
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
                "IMPACT_filename": {
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
                    'checksum': 'sha1$6708055e2cadbbb9d409d828ec0d6c37d340274f',
                    'size': 70,
                    'path': os.path.join(output_dir, 'output.maf')
                    }
                }
            self.assertDictEqual(output_json, expected_output)

            comments, mutations = load_mutations(output_json['IMPACT_col_added_file']['path'])

            expected_comments = ['# comment 1', '# comment 2']
            self.assertEqual(comments, expected_comments)

            expected_mutations = [
                {'Hugo_Symbol': 'SUFU', 'is_in_impact': 'True'},
                {'Hugo_Symbol': 'GOT1', 'is_in_impact': 'False'}
                ]
            self.assertEqual(mutations, expected_mutations)

    def test_add_impact_1(self):
        """
        Test that a maf file with is_in_IMPACT column comes out as expected
        """
        self.maxDiff = None
        input_maf = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf")
        impact_file = os.path.join(IMPACT_FILE)

        with TemporaryDirectory() as tmpdir:
            # output_dir = os.path.join(tmpdir, "output")
            input_json = {
                "input_file": {
                      "class": "File",
                      "path": input_maf
                    },
                "output_filename":  'output.maf',
                "IMPACT_filename": {
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
                    'checksum': 'sha1$d94362b7f33c11c9b64bcc68ee356e933da56c1d',
                    'size': 15323931,
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

            # check that its got an extra column in the output
            self.assertTrue(len(input_mutations[1])+1==len(output_mutations[1]))

if __name__ == "__main__":
    unittest.main()
