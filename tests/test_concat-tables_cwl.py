#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the concat-tables.cwl
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import CWLFile, PlutoTestCase
sys.path.pop(0)

class TestConcatTables(PlutoTestCase):
    cwl_file = CWLFile('concat-tables.cwl')
    def test_concat_two_tables(self):
        """
        Test that two files are concatenated correctly
        """
        # make a dummy file with some lines
        input_lines1 = ["HEADER1", "foo1", "bar1"]
        input_file1 = os.path.join(self.tmpdir, "input1.txt")
        with open(input_file1, "w") as fout:
            for line in input_lines1:
                fout.write(line + '\n')

        input_lines2 = ["HEADER2", "foo2", "bar2"]
        input_file2 = os.path.join(self.tmpdir, "input2.txt")
        with open(input_file2, "w") as fout:
            for line in input_lines2:
                fout.write(line + '\n')

        self.input = {
            "input_files": [{
                  "class": "File",
                  "path": input_file1
                },
                {
                  "class": "File",
                  "path": input_file2
                }
                ],
            "output_filename": "output.txt"
            }

        output_json, output_dir = self.run_cwl()

        # check the contents of the concatenated file; should be the same as the input
        output_file = output_json['output_file']['path']
        with open(output_file) as fin:
            output_lines = [ line.strip() for line in fin ]

        expected_lines = ['HEADER1\tHEADER2', 'foo1\tNA', 'bar1\tNA', 'NA\tfoo2', 'NA\tbar2']
        self.assertEqual(output_lines, expected_lines)

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir, 'output.txt'),
                'basename': 'output.txt',
                'class': 'File',
                'checksum': 'sha1$d92a4e707cb5dad2ec557edfe976680dfffc5f3f',
                'size': 53,
                'path': os.path.join(output_dir, 'output.txt')
                }
            }
        self.assertDictEqual(output_json, expected_output)

    def test_concat_one_tables(self):
        """
        Test that one file is returned correctly from the script
        """
        # make a dummy file with some lines
        input_lines1 = ["HEADER1", "foo1", "bar1"]
        input_file1 = os.path.join(self.tmpdir, "input1.txt")
        with open(input_file1, "w") as fout:
            for line in input_lines1:
                fout.write(line + '\n')

        self.input = {
            "input_files": [{
                  "class": "File",
                  "path": input_file1
                },
                ],
            "output_filename": "output.txt"
            }

        output_json, output_dir = self.run_cwl()

        # check the contents of the concatenated file; should be the same as the input
        output_file = output_json['output_file']['path']
        with open(output_file) as fin:
            output_lines = [ line.strip() for line in fin ]

        expected_lines = ['HEADER1', 'foo1', 'bar1']
        self.assertEqual(output_lines, expected_lines)

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir, 'output.txt'),
                'basename': 'output.txt',
                'class': 'File',
                'checksum': 'sha1$2274c54c24a98e8235e34d78b700d04cb95f48dd',
                'size': 21,
                'path': os.path.join(output_dir, 'output.txt')
                }
            }
        self.assertDictEqual(output_json, expected_output)


if __name__ == "__main__":
    unittest.main()
