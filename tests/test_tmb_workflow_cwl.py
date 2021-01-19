#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the TMB analysis workflow cwl which uses multiple input samples
"""
import os
import unittest

# relative imports, from CLI and from parent project
if __name__ != "__main__":
    from .tools import TmpDirTestCase, load_mutations, run_cwl, write_table, dicts2lines
    from .settings import CWL_DIR

if __name__ == "__main__":
    from tools import TmpDirTestCase, load_mutations, run_cwl, write_table, dicts2lines
    from settings import CWL_DIR

cwl_file = os.path.join(CWL_DIR, 'tmb_workflow.cwl')

class TestTmbWorkflow(TmpDirTestCase):
    def setUp(self):
        # initialize the tmpdir
        super().setUp()
        self.data_clinical_lines = [
        ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE'],
        ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE'],
        ['#STRING', 'STRING', 'NUMBER',],
        ['#1', '1', '1'],
        ['SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE'],
        ['Sample1-T', 'Patient1', '108'],
        ['Sample1-N', 'Patient2', '58'],
        ['Sample2-T', 'Patient3', '502'],
        ['Sample2-N', 'Patient4', '56'],
        ]
        self.maf_comments = [
        ['# comment 1'],
        ['# comment 2']
        ]
        self.maf_row1 = {
        't_af': '0.50',
        't_depth': '550',
        'Hugo_Symbol': 'EGFR',
        'Start_Position': '1',
        'Consequence': 'synonymous_variant' # exclude due to synonymous_variant
        }
        self.maf_row1_2 = {
        't_af': '0.51',
        't_depth': '551',
        'Hugo_Symbol': 'EGFR',
        'Start_Position': '1',
        'Consequence': 'synonymous_variant' # exclude due to synonymous_variant
        }
        self.maf_row2 = {
        't_af': '0.50',
        't_depth': '550',
        'Hugo_Symbol': 'EGFR',
        'Start_Position': '1',
        'Consequence': 'splice_region_variant,synonymous_variant' # exclude due to synonymous_variant
        }
        self.maf_row3 = { # this one should pass filter
        't_af': '0.50',
        't_depth': '550',
        'Hugo_Symbol': 'EGFR',
        'Start_Position': '1',
        'Consequence': 'missense_variant'
        }
        self.maf_row3_2 = { # this one should pass filter
        't_af': '0.52',
        't_depth': '552',
        'Hugo_Symbol': 'EGFR',
        'Start_Position': '1',
        'Consequence': 'missense_variant'
        }
        self.maf_row3_3 = {
        't_af': '0.01', # exclude due to low AF
        't_depth': '552',
        'Hugo_Symbol': 'EGFR',
        'Start_Position': '1',
        'Consequence': 'missense_variant'
        }
        self.maf_row4 = {
        't_af': '0.01', # exclude due to low AF
        't_depth': '550',
        'Hugo_Symbol': 'EGFR',
        'Start_Position': '1',
        'Consequence': 'missense_variant'
        }
        self.maf_row5 = {
        't_af': '0.51',
        't_depth': '90', # exclude due to low coverage
        'Hugo_Symbol': 'EGFR',
        'Start_Position': '1',
        'Consequence': 'missense_variant'
        }
        self.maf_row6 = { # this one should pass filter
        't_af': '0.45',
        't_depth': '590',
        'Hugo_Symbol': 'EGFR',
        'Start_Position': '1',
        'Consequence': 'splice_region_variant'
        }
        self.maf_row7 = { # this one should pass filter
        't_af': '0.45',
        't_depth': '590',
        'Hugo_Symbol': 'TERT',
        'Start_Position': '1295340', # good value; is_TERT_promoter = True
        'Consequence': 'splice_region_variant'
        }
        self.maf_row8 = { # this should pass filter
        't_af': '0.45',
        't_depth': '590',
        'Hugo_Symbol': 'TERT',
        'Start_Position': '1295339', # good value; is_TERT_promoter = True
        'Consequence': 'splice_region_variant'
        }
        self.maf_row9 = { # this should pass filter
        't_af': '0.45',
        't_depth': '590',
        'Hugo_Symbol': 'TERT',
        'Start_Position': '1295341', # bad value; is_TERT_promoter = False
        'Consequence': 'splice_region_variant' # include anyway because its not synonymous_variant
        }
        self.maf_row10 ={ # this should pass filter
        't_af': '0.45',
        't_depth': '590',
        'Hugo_Symbol': 'TERT',
        'Start_Position': '1295339', # good value; is_TERT_promoter = True
        'Consequence': 'synonymous_variant' # include even though its synonymous_variant
        }
        self.maf_rows1 = [
            self.maf_row1,
            self.maf_row2,
            self.maf_row3,
            self.maf_row4,
            self.maf_row5,
            self.maf_row6,
            self.maf_row7,
            self.maf_row8,
            self.maf_row9,
            self.maf_row10
        ]
        self.maf_rows2 = [
            self.maf_row1_2,
            self.maf_row2,
            self.maf_row3_2,
            self.maf_row3_3,
            self.maf_row4,
            self.maf_row5,
            self.maf_row6,
            self.maf_row8,
            self.maf_row9,
            self.maf_row10
        ]
        self.maf_lines1 = dicts2lines(dict_list = self.maf_rows1, comment_list = self.maf_comments)
        self.maf_lines2 = dicts2lines(dict_list = self.maf_rows2, comment_list = self.maf_comments)
        self.maf1 = write_table(self.tmpdir, filename = "input1.maf", lines = self.maf_lines1)
        self.maf2 = write_table(self.tmpdir, filename = "input2.maf", lines = self.maf_lines2)
        self.data_clinical_file = write_table(self.tmpdir, filename = "data_clinical_sample.txt", lines = self.data_clinical_lines)

    def test_tmb_workflow1(self):
        """
        Test case for running the TMB workflow on multiple files
        """
        self.maxDiff = None
        input_json = {
            "data_clinical_file": {
                  "class": "File",
                  "path": self.data_clinical_file
                },
            "assay_coverage":  '1000',
            "pairs": [
                {
                    "pair_maf": {
                        "path": self.maf1,
                        "class": "File"
                    },
                    "pair_id": "Sample1-T.Sample1-N",
                    "tumor_id": "Sample1-T",
                    "normal_id": "Sample1-N"
                },
                {
                    "pair_maf": {
                        "path": self.maf2,
                        "class": "File"
                    },
                    "pair_id": "Sample2-T.Sample2-N",
                    "tumor_id": "Sample2-T",
                    "normal_id": "Sample2-N"
                }
                ]
            }
        output_json, output_dir = run_cwl(
            testcase = self,
            tmpdir = self.tmpdir,
            input_json = input_json,
            cwl_file = cwl_file,
            print_command = False,
            )

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'data_clinical_sample.txt'),
                'basename': 'data_clinical_sample.txt',
                'class': 'File',
                'checksum': 'sha1$6933552a71f7d08bb35f685ce05565399a731c8d',
                'size': 307,
                'path':  os.path.join(output_dir,'data_clinical_sample.txt')
                }
            }
        self.assertDictEqual(output_json, expected_output)

        output_file = expected_output['output_file']['path']
        with open(output_file) as fin:
            lines = [ l.strip().split() for l in fin ]

        expected_lines = [
            ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_TMB_SCORE'],
            ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_TMB_SCORE'],
            ['#STRING', 'STRING', 'NUMBER', 'NUMBER'],
            ['#1', '1', '1', '1'],
            ['SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_TMB_SCORE'],
            ['Sample1-T', 'Patient1', '108', '0.000000006'],
            ['Sample1-N', 'Patient2', '58'], # tailing empty value gets stripped off
            ['Sample2-T', 'Patient3', '502', '0.000000005'],
            ['Sample2-N', 'Patient4', '56'] # tailing empty value gets stripped off
            ]
        self.assertEqual(lines, expected_lines)

if __name__ == "__main__":
    unittest.main()
