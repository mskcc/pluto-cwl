#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for merging tables
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import TmpDirTestCase, run_cwl, write_table, CWLFile

sys.path.pop(0)

cwl_file = CWLFile('merge-tables.cwl')

class TestMergeTables(TmpDirTestCase):
    def test_merge_tables1(self):
        lines1 = [
        ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE'],
        ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE'],
        ['#STRING', 'STRING', 'NUMBER'],
        ['#1', '1', '1'],
        ['SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE'],
        ['Sample1-T', 'Patient1', '108'],
        ['Sample1-N', 'Patient2', '58'],
        ['Sample2-T', 'Patient3', '502'],
        ['Sample2-N', 'Patient4', '56'],
        ]

        lines2 = [
        ['SampleID', 'CMO_TMB_SCORE'],
        ['Sample1-T', '100'],
        ['Sample2-T', '200'],
        ]
        data_clinical_file = write_table(self.tmpdir, filename = "data_clinical_sample.txt", lines = lines1)
        tmb_file = write_table(self.tmpdir, filename = "tmb.tsv", lines = lines2)
        input_json = {
            'table1': {
                "class": "File",
                "path": data_clinical_file
            },
            'table2': {
                "class": "File",
                "path": tmb_file
            },
            'key1': 'SAMPLE_ID',
            'key2': 'SampleID',
            'output_filename': 'output.tsv',
            'cBioPortal': True
            #  NOTE: need to use cBioPortal flag to get correct output headers
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
                'location': 'file://' + os.path.join(output_dir,'output.tsv'),
                'basename': 'output.tsv',
                'class': 'File',
                'checksum': 'sha1$03d7120a05eeef1a4a7349c9b35069566415e786',
                'size': 255,
                'path':  os.path.join(output_dir,'output.tsv')
                }
            }
        self.assertDictEqual(output_json, expected_output)

        output_file = expected_output['output_file']['path']
        with open(output_file) as fin:
            lines = [ l.strip().split() for l in fin ]

        #  NOTE: need to use cBioPortal flag to get correct output headers
        expected_lines = [
            ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_TMB_SCORE'],
            ['#SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_TMB_SCORE'],
            ['#STRING', 'STRING', 'NUMBER', 'NUMBER'],
            ['#1', '1', '1', '1'],
            ['SAMPLE_ID', 'PATIENT_ID', 'SAMPLE_COVERAGE', 'CMO_TMB_SCORE'],
            ['Sample1-T', 'Patient1', '108', '100'],
            ['Sample1-N', 'Patient2', '58', 'NA'],
            ['Sample2-T', 'Patient3', '502', '200'],
            ['Sample2-N', 'Patient4', '56', 'NA']
            ]
        self.assertEqual(lines, expected_lines)

if __name__ == "__main__":
    unittest.main()
