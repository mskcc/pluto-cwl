#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the snp-pileup-wrapper.cwl file
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import CWLFile, PlutoTestCase
from pluto.serializer import OFile
sys.path.pop(0)

class TestConcatWithCommentsCWL(PlutoTestCase):
    cwl_file = CWLFile('concat_with_comments.cwl')

    def test_concat_0(self):
        """
        Test concat when no comments are present in the original file
        """
        # make a dummy file with some lines
        input_lines = ["HEADER", "foo", "bar", "baz"]
        input_file = os.path.join(self.tmpdir, "input.txt")
        with open(input_file, "w") as fout:
            for line in input_lines:
                fout.write(line + '\n')

        self.input = {
            "input_files": [{
                  "class": "File",
                  "path": input_file
                }],
            "comment_label": "comment_label",
            "comment_value": "comment_value"
            }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': OFile(name='output.txt', size=49, hash='7cef8f6de47289a55de99de77563beb3fa371deb', dir = output_dir)
            }
        self.assertCWLDictEqual(output_json, expected_output)

        # check the contents of the concatenated file; should be the same as the input
        self.assertFileLinesEqual(
            expected_output['output_file']['path'],
            ['#comment_label: comment_value', "HEADER", 'foo', 'bar', 'baz'])


    def test_concat1(self):
        """
        Test concat when original file has a comment line
        """
        # make a dummy file with some lines
        input_lines = ["# comment here", "HEADER", "foo", "bar", "baz"]
        input_file = os.path.join(self.tmpdir, "input.txt")
        with open(input_file, "w") as fout:
            for line in input_lines:
                fout.write(line + '\n')

        self.input = {
            "input_files": [{
                  "class": "File",
                  "path": input_file
                }],
            "comment_label": "comment_label",
            "comment_value": "comment_value"
            }
        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': OFile(name='output.txt', size=64, hash='14ee1247f314dba1e3c28aa8aec9ff7b137a1f41', dir = output_dir)
            }
        self.assertCWLDictEqual(output_json, expected_output)

        # check the contents of the concatenated file; should be the same as the input
        self.assertFileLinesEqual(
            expected_output['output_file']['path'],
            ['# comment here', '#comment_label: comment_value', "HEADER", 'foo', 'bar', 'baz'])


    def test_concat2(self):
        """
        Test concat when multiple files have comments
        """
        # make a dummy file with some lines
        input_lines1 = ["# comment 1 here", "HEADER", "foo1", "bar1"]
        input_file1 = os.path.join(self.tmpdir, "input1.txt")
        with open(input_file1, "w") as fout:
            for line in input_lines1:
                fout.write(line + '\n')

        input_lines2 = ["# comment 2 here", "HEADER", "foo2", "bar2"]
        input_file2 = os.path.join(self.tmpdir, "input2.txt")
        with open(input_file2, "w") as fout:
            for line in input_lines2:
                fout.write(line + '\n')

        self.input = {
            "input_files": [
                {
                  "class": "File",
                  "path": input_file1
                },
                {
                  "class": "File",
                  "path": input_file2
                },
                ],
            "comment_label": "comment_label",
            "comment_value": "comment_value"
            }
        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': OFile(name='output.txt', size=91, hash='5dbce16f9bfef135d6b8288b16350351a33998f3', dir = output_dir)
            }
        self.assertCWLDictEqual(output_json, expected_output)

        self.assertFileLinesEqual(
            expected_output['output_file']['path'],
            [
                '# comment 1 here',
                '# comment 2 here',
                '#comment_label: comment_value',
                "HEADER",
                'foo1',
                'bar1',
                'foo2',
                'bar2'
                ])

if __name__ == "__main__":
    unittest.main()
