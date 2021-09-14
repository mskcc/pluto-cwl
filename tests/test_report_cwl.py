#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for making a report
"""
import os
import sys
import unittest

PARENT_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile
from pluto.serializer import OFile
sys.path.pop(0)

class TestMergeTables(PlutoTestCase):
    cwl_file = CWLFile('report.cwl')

    def test_demo_report(self):
        """
        Test case for making an HTML report
        """
        mut1 = {
        "Hugo_Symbol" : "FGF3",
        "Amino_Acid_Change": "E116K",
        "Tumor_Sample_Barcode": "Sample1",
        "Matched_Norm_Sample_Barcode": "Sample2",
        "t_depth" : "0",
        "t_alt_count" : "0",
        "n_alt_count": "0",
        "n_depth": "0"
        }

        sample_lines = [
        {
        'SAMPLE_ID': 'Sample1',
        'PATIENT_ID': 'Patient1',
        'IGO_ID': 'IGO_1',
        'COLLAB_ID': 'Collab1',
        'ONCOTREE_CODE': 'MEL',
        'SAMPLE_COVERAGE': '100',
        },
        {
        'SAMPLE_ID': 'Sample2',
        'PATIENT_ID': 'Patient1',
        'IGO_ID': 'IGO_2',
        'COLLAB_ID': 'Collab2',
        'ONCOTREE_CODE': 'MEL',
        'SAMPLE_COVERAGE': '100',
        },
        ]

        patient_lines = [
        ['#PATIENT_ID', 'SEX'],
        ['#PATIENT_ID', 'SEX'],
        ['#STRING', 'STRING'],
        ['#1', '1'],
        ['PATIENT_ID', 'SEX'],
        ['p_C_00001', 'F']
        ]

        mut_lines = self.dicts2lines(dict_list = [mut1], comment_list = [])
        mut_file = self.write_table(self.tmpdir, filename = "input.maf", lines = mut_lines)

        sample_lines = self.dicts2lines(dict_list = sample_lines, comment_list = [])
        sample_file = self.write_table(self.tmpdir, filename = "samples.txt", lines = sample_lines)

        patient_file = self.write_table(self.tmpdir, filename = "patients.txt", lines = patient_lines)
        output_file = os.path.join(self.tmpdir, "report.html")

        self.input = {
            'mutation_file': { "class": "File", "path": mut_file },
            'samples_file': { "class": "File", "path": sample_file },
            'patients_file': { "class": "File", "path": patient_file },
            'output_filename': 'report.html'
        }

        output_json, output_dir = self.run_cwl()

        # HTML output files have embedded timestamp that change size and hash
        output_json['output_file'].pop('checksum')
        output_json['output_file'].pop('size')

        expected_output = {
            'output_file': OFile(name = 'report.html', dir = output_dir)
            }

        self.assertDictEqual(output_json, expected_output)

if __name__ == "__main__":
    unittest.main()