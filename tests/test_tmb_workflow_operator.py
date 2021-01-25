#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the TMB analysis workflow Operator
"""
import os
import sys
import unittest
import json

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import TmpDirTestCase, TableReader, load_mutations, run_cwl, write_table, dicts2lines, CWLFile
from operators.tmb_workflow import TMBWorkflow
sys.path.pop(0)

class TestTmbWorkflowOperator(TmpDirTestCase):
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
        ['Sample6-T', 'Patient4', '57'],
        ['Sample7-N', 'Patient4', '58'],
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
            self.maf_row1, # exclude due to synonymous_variant
            self.maf_row2, # exclude due to synonymous_variant
            self.maf_row3, # this one should pass filter
            self.maf_row4, # exclude due to low AF
            self.maf_row5, # exclude due to low coverage
            self.maf_row6, # this one should pass filter
            self.maf_row7,  # this one should pass filter
            self.maf_row8, # this should pass filter
            self.maf_row9, # this should pass filter
            self.maf_row10 # this should pass filter
        ]
        # PASS: 6
        # FAIL: 4
        self.maf_rows2 = [
            self.maf_row1_2, # exclude due to synonymous_variant
            self.maf_row2, # exclude due to synonymous_variant
            self.maf_row3_2, # this one should pass filter
            self.maf_row3_3, # exclude due to low AF
            self.maf_row4, # exclude due to low AF
            self.maf_row5, # exclude due to low coverage
            self.maf_row6, # this one should pass filter
            self.maf_row8, # this should pass filter
            self.maf_row9, # this should pass filter
            self.maf_row10 # this should pass filter
        ]
        # PASS: 5
        # FAIL: 5
        self.maf_lines1 = dicts2lines(dict_list = self.maf_rows1, comment_list = self.maf_comments)
        self.maf_lines2 = dicts2lines(dict_list = self.maf_rows2, comment_list = self.maf_comments)
        self.maf1 = write_table(self.tmpdir, filename = "input1.maf", lines = self.maf_lines1)
        self.maf2 = write_table(self.tmpdir, filename = "input2.maf", lines = self.maf_lines2)
        self.data_clinical_file = write_table(self.tmpdir, filename = "data_clinical_sample.txt", lines = self.data_clinical_lines)

        self.pairs_dicts = [
            {
                'tumor_id': 'Sample1-T',
                'normal_id': 'Sample1-N',
                'pair_id': 'Sample1-T.Sample1-N',
                'pair_maf': self.maf1,
                'snp_pileup': ''
            },
            {
                'tumor_id': 'Sample2-T',
                'normal_id': 'Sample2-N',
                'pair_id': 'Sample2-T.Sample2-N',
                'pair_maf': self.maf2,
                'snp_pileup': ''
            }
        ]
        self.pairs_lines = dicts2lines(dict_list = self.pairs_dicts, comment_list = [])
        self.pairs_file = write_table(self.tmpdir, filename = "pairs.tsv", lines = self.pairs_lines)

    def test_tmb_workflow_operator1(self):
        """
        Test case for running the TMB operator with
        """
        self.maxDiff = None
        operator = TMBWorkflow(data_clinical_file = self.data_clinical_file,
            assay_coverage = '1000',
            pairs_file = self.pairs_file,
            dir = self.tmpdir,
            verbose = False)
        output_json, output_dir, output_json_file = operator.run()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'data_clinical_sample.txt'),
                'basename': 'data_clinical_sample.txt',
                'class': 'File',
                'checksum': 'sha1$e7975b7d9f34202750c4645d02f2f98796e22c74',
                'size': 363,
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
            ['Sample1-N', 'Patient2', '58', 'NA'],
            ['Sample2-T', 'Patient3', '502', '0.000000005'],
            ['Sample2-N', 'Patient4', '56', 'NA'],
            ['Sample6-T', 'Patient4', '57', 'NA'],
            ['Sample7-N', 'Patient4', '58', 'NA']
            ]
        self.assertEqual(lines, expected_lines)

if __name__ == "__main__":
    unittest.main()
