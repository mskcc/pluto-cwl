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
from pluto.tools import PlutoTestCase, CWLFile
sys.path.pop(0)

class TestAddMSIStatus(PlutoTestCase):
    cwl_file = CWLFile('add_msi_status.cwl')

    def test_add_msi_status(self):
        """
        Test case for adding the MSI Status label to a table
        """
        lines1 = [
            ['Total_Number_of_Sites', 'Number_of_Somatic_Sites', 'MSI_SCORE', 'SAMPLE_ID'],
            ['123',                   '987',                     '11',     'Sample1-T'],
            ['456',                   '654',                     '2',      'Sample2-T'],
            ['789',                   '321',                     '5',      'Sample3-T']
        ]

        tmb_file = self.write_table(self.tmpdir, filename = "msi.tsv", lines = lines1)
        self.input = {
            'input_filename': {
                "class": "File",
                "path": tmb_file
            },
            'output_filename': 'output.tsv',
            'header': 'MSI_STATUS'
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'output.tsv'),
                'basename': 'output.tsv',
                'class': 'File',
                'checksum': 'sha1$ea9eb5b227f25bd03a50fbcec282259e89d176ac',
                'size': 168,
                'path':  os.path.join(output_dir,'output.tsv')
                }
            }
        self.assertDictEqual(output_json, expected_output)

        output_file = expected_output['output_file']['path']

        lines = self.read_table(output_file)

        expected_lines = [
            ['Total_Number_of_Sites',   'Number_of_Somatic_Sites', 'MSI_SCORE',       'SAMPLE_ID',       'MSI_STATUS'],
            ['123',                     '987',                     '11',              'Sample1-T',       'Instable'],
            ['456',                     '654',                     '2',               'Sample2-T',       'Stable'],
            ['789',                     '321',                     '5',               'Sample3-T',       'Indeterminate']
            ]
        self.assertEqual(lines, expected_lines)



if __name__ == "__main__":
    unittest.main()
