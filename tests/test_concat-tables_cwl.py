#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the concat-tables.cwl
"""
import os
import sys
import unittest



from pluto import (
    CWLFile, 
    PlutoTestCase,
)


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
        output_path = os.path.join(output_dir, 'output.txt')

        # check the contents of the concatenated file; should be the same as the input
        output_file = output_path
        with open(output_file) as fin:
            output_lines = [ line.strip() for line in fin ]

        expected_lines = ['HEADER1\tHEADER2', 'foo1\tNA', 'bar1\tNA', 'NA\tfoo2', 'NA\tbar2']
        self.assertEqual(output_lines, expected_lines)

        expected_output = {
            'output_file': {
                'location': 'file://' + output_path,
                'basename': 'output.txt',
                'class': 'File',
                'checksum': 'sha1$1d1d9666cf862439ff6e27c3735325dab73c5830',
                'size': 48,
                'path': output_path
                }
            }
        self.assertCWLDictEqual(output_json, expected_output)

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
        output_path = os.path.join(output_dir, 'output.txt')

        # check the contents of the concatenated file; should be the same as the input
        output_file = output_path
        with open(output_file) as fin:
            output_lines = [ line.strip() for line in fin ]

        expected_lines = ['HEADER1', 'foo1', 'bar1']
        self.assertEqual(output_lines, expected_lines)

        expected_output = {
            'output_file': {
                'location': 'file://' + output_path,
                'basename': 'output.txt',
                'class': 'File',
                'checksum': 'sha1$83839e4ffe066577c770d8fd0a530530899a4751',
                'size': 18,
                'path': output_path
                }
            }
        self.assertCWLDictEqual(output_json, expected_output)




