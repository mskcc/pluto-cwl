#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for merging tables
"""
import os
import sys
import unittest

PARENT_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, PARENT_DIR)
from pluto import (
    PlutoTestCase, 
    CWLFile
)
sys.path.pop(0)

class TestMergeTables(PlutoTestCase):
    cwl_file = CWLFile('merge-tables.cwl')

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
        data_clinical_file = self.write_table(self.tmpdir, filename = "data_clinical_sample.txt", lines = lines1)
        tmb_file = self.write_table(self.tmpdir, filename = "tmb.tsv", lines = lines2)
        self.input = {
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

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'output.tsv'),
                'basename': 'output.tsv',
                'class': 'File',
                'checksum': 'sha1$d789b7d7e9c6e60afc37c907979e52b8e9b48841',
                'size': 297,
                'path':  os.path.join(output_dir,'output.tsv')
                }
            }
        self.assertCWLDictEqual(output_json, expected_output)

        output_file = os.path.join(output_dir,'output.tsv')

        lines = self.read_table(output_file)

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

    def test_merge_tables2(self):
        lines1 = [
        ["SampleID", "Foo"],
        ["A", '1'],
        ["B", '2']
        ]
        lines2 = [
        ["SampleID", "Bar"],
        ["A", '100'],
        ["B", '200']
        ]
        file1 = self.write_table(self.tmpdir, filename = "file1.txt", lines = lines1)
        file2 = self.write_table(self.tmpdir, filename = "file2.tsv", lines = lines2)
        self.input = {
            'table1': {
                "class": "File",
                "path": file1
            },
            'table2': {
                "class": "File",
                "path": file2
            },
            'key1': 'SampleID',
            'key2': 'SampleID',
            'output_filename': 'output.tsv'
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'output.tsv'),
                'basename': 'output.tsv',
                'class': 'File',
                'checksum': 'sha1$5a770c640ff64a5c84f13b77b58c0b03fb18bd97',
                'size': 36,
                'path':  os.path.join(output_dir,'output.tsv')
                }
            }

        self.assertCWLDictEqual(output_json, expected_output)

        output_file = os.path.join(output_dir,'output.tsv')

        lines = self.read_table(output_file)

        expected_lines = [
            ['SampleID', 'Foo', 'Bar'],
            ['A', '1', '100'],
            ['B', '2', '200']
            ]
        self.assertEqual(lines, expected_lines)

if __name__ == "__main__":
    unittest.main()
