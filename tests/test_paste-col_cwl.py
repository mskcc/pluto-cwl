#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the paste-col.cwl
"""
import os
import sys
import unittest



from pluto import (
    PlutoTestCase, 
    CWLFile,
    OFile
)


class TestPasteCol(PlutoTestCase):
    cwl_file = CWLFile('paste-col.cwl')

    def test_paste_col_1(self):
        """
        """
        # make a dummy file with some lines
        input_lines = ["HEADER1", "foo1", "bar1"]
        input_file = os.path.join(self.tmpdir, "input.txt")
        with open(input_file, "w") as fout:
            for line in input_lines:
                fout.write(line + '\n')

        self.input = {
            "input_file": {
                  "class": "File",
                  "path": input_file
                },
            "output_filename": "output.txt",
            "header": "HEADER2",
            "value": "foo2"
            }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': OFile(name = "output.txt", size = 36, hash = "34753fd98b2355d54740f3fdfc6490262c15dd59", dir = output_dir)
            }

        self.assertCWLDictEqual(output_json, expected_output)

        output_file = expected_output['output_file']['path']
        with open(output_file) as fin:
            output_lines = [ line.strip() for line in fin ]

        expected_lines = ['HEADER1\tHEADER2', 'foo1\tfoo2', 'bar1\tfoo2']
        self.assertEqual(output_lines, expected_lines)




