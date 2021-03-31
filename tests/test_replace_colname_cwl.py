#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for CWL to change header column name
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile, TableReader
from pluto.settings import DATA_SETS
sys.path.pop(0)

class TestReplaceColname(PlutoTestCase):
    cwl_file = CWLFile('replace_colname.cwl')

    def setUp(self):
        # initialize the tmpdir
        super().setUp()
        lines = [
        ["A", "%", "C"],
        ["1", "2", "3"],
        ["foo", "bar", "baz"]
        ]
        self.table = self.write_table(self.tmpdir, filename = "table.tsv", lines = lines)

    def test_change_colnames(self):
        self.maxDiff = None
        self.input = {
            "old_name": "%",
            "new_name": "pcnt",
            "input_file": {"class": "File", "path": self.table}
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'output.tsv'),
                'basename': 'output.tsv',
                'class': 'File',
                'checksum': 'sha1$dfe58e029a2c3262dbdf4270b0b5af6bf420589e',
                'size': 27,
                'path':  os.path.join(output_dir,'output.tsv')
                }
            }
        self.assertDictEqual(output_json, expected_output)

        output_file = expected_output['output_file']['path']

        lines = self.read_table(output_file)

        expected_lines = [
        ['A', 'pcnt', 'C'],
        ['1', '2', '3'],
        ['foo', 'bar', 'baz']
        ]
        self.assertEqual(lines, expected_lines)

if __name__ == "__main__":
    unittest.main()
