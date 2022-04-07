#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile

sys.path.pop(0)

class TestAddHeader(PlutoTestCase):
    cwl_file = CWLFile('add_header.cwl')

    def test_add_header(self):
        """
        Test case for adding a header to a file
        """
        self.maxDiff = None

        input_file = os.path.join(self.tmpdir, "input.txt")
        with open(input_file, "w") as f:
            f.write("foo")
        header_str = "HEADER"

        self.input = {
            "input_file": {
                  "class": "File",
                  "path": input_file
                },
            "header_str":  header_str,
            }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'output.txt'),
                'basename': 'output.txt',
                'class': 'File',
                'checksum': 'sha1$01838a0977d542fb12680e271393e1d4baaefa8f',
                'size': 10,
                'path':  os.path.join(output_dir,'output.txt')
                }
            }
        self.assertCWLDictEqual(output_json, expected_output)

        output_file = expected_output['output_file']['path']
        with open(output_file) as f:
            lines = [ l.strip() for l in f ]
        expected_lines = ['HEADER', 'foo']
        self.assertEqual(lines, expected_lines)

    def test_add_header_empty_file(self):
        """
        Test case for adding a header to an empty file should return only the header
        """
        input_file = os.path.join(self.tmpdir, "input.txt")
        with open(input_file, "w") as f:
            pass

        header_str = "HEADER"

        self.input = {
            "input_file": {
                  "class": "File",
                  "path": input_file
                },
            "header_str":  header_str,
            }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': {
                'location': 'file://' + os.path.join(output_dir,'output.txt'),
                'basename': 'output.txt',
                'class': 'File',
                'checksum': 'sha1$b4cf58442d6321c81db6bab562806e14bf54bf72',
                'size': 7,
                'path':  os.path.join(output_dir,'output.txt')
                }
            }
        self.assertCWLDictEqual(output_json, expected_output)

        output_file = expected_output['output_file']['path']
        with open(output_file) as f:
            lines = [ l.strip() for l in f ]
        expected_lines = ['HEADER']
        self.assertEqual(lines, expected_lines)

if __name__ == "__main__":
    unittest.main()
