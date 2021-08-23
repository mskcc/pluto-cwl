#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for merging CNA tables via the full outer join CWL
"""
import os
import sys
import unittest

PARENT_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile
sys.path.pop(0)

class TestMergeTables(PlutoTestCase):
    cwl_file = CWLFile('full-outer-join.cwl')

    def test_merge_tables1(self):
        lines1 = [
        ['Hugo_Symbol', 'Sample1', 'Sample2'],
        ["TAP1", "0", "0"],
        ["ERRFI1", "0", "0"],
        ["STK19", "", "0"],
        ]

        lines2 = [
        ['Hugo_Symbol', 'Sample3', 'Sample4'],
        ["ERRFI1", "0", "0"],
        ["STK19", "-2", "0"],
        ["STK11", "0", ""],
        ]
        cna_file1 = self.write_table(self.tmpdir, filename = "cna1.txt", lines = lines1)
        cna_file2 = self.write_table(self.tmpdir, filename = "cna2.txt", lines = lines2)
        self.input = {
            'table1': { "class": "File", "path": cna_file1 },
            'table2': { "class": "File", "path": cna_file2 },
            'join_key': 'Hugo_Symbol',
            'output_filename': 'output.tsv'
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'output.tsv'),
                'basename': 'output.tsv',
                'class': 'File',
                'checksum': 'sha1$1548027e2d346ccad32414e0fe35d215d519eda1',
                'size': 107,
                'path':  os.path.join(output_dir,'output.tsv')
                }
            }
        self.assertDictEqual(output_json, expected_output)

        output_file = expected_output['output_file']['path']

        lines = self.read_table(output_file)

        expected_lines = [
            ['Hugo_Symbol', 'Sample1', 'Sample2', 'Sample3', 'Sample4'],
            ['ERRFI1', '0', '0', '0', '0'],
            ['STK19', 'NA', '0', '-2', '0'],
            ['TAP1', '0', '0', 'NA', 'NA'],
            ['STK11', 'NA', 'NA', '0', 'NA']
            ]
        self.assertEqual(lines, expected_lines)


if __name__ == "__main__":
    unittest.main()
