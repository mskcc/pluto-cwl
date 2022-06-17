#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the replace.cwl file
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile
from pluto.serializer import OFile
sys.path.pop(0)

class TestReplace(PlutoTestCase):
    cwl_file = CWLFile('replace.cwl')

    def test_replace1(self):
        """
        Test that strings get replaced
        """
        # make a dummy file with some lines
        input_lines = ["HEADER", "foo", "ILLOGICAL", "baz"]
        input_file = os.path.join(self.tmpdir, "input.txt")
        with open(input_file, "w") as fout:
            for line in input_lines:
                fout.write(line + '\n')

        self.input = {
            "input_file":{
                  "class": "File",
                  "path": input_file
                },
            }

        output_json, output_dir = self.run_cwl()

        # check the contents of the concatenated file; should be the same as the input
        output_file = os.path.join(output_dir, 'output.txt')
        with open(output_file) as fin:
            output_lines = [ line.strip() for line in fin ]

        expected_lines = ["HEADER", "foo", "NA", "baz"]
        self.assertEqual(output_lines, expected_lines)

        expected_output = {
            'output_file': OFile(name = "output.txt", size = 18, hash = "62255c8ee13b8ba6e01c7e17262a8ba1f174e5cb", dir = output_dir)
            }
        self.assertCWLDictEqual(output_json, expected_output)


if __name__ == "__main__":
    unittest.main()
