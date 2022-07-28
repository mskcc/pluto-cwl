#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unit tests for reduce_sig_figs.cwl
"""
import os
import sys
import unittest
import csv
from collections import OrderedDict

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import CWLFile, PlutoTestCase
from pluto.serializer import OFile
sys.path.pop(0)

class TestReduceSigFigs(PlutoTestCase):
    cwl_file = CWLFile('reduce_sig_figs.cwl')

    def test_reduce_sig_figs(self):
        """
        Test that significant figures are reduced correctly
        """
        input_lines = ["seg.mean", "3.141592", "2.718281828"]
        input_file = os.path.join(self.tmpdir, "input.txt")
        with open(input_file, "w") as fout:
            for line in input_lines:
                fout.write(line + '\n')

        self.input = {
            "input_file": {
                  "class": "File",
                  "path": input_file
                }
            }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'output_file': OFile(name = "output.txt", size = 26, hash = "d9f5ec4a9aa27a69ee64edb97eb10d6db65c7ad7", dir = output_dir)
        }

        self.assertCWLDictEqual(output_json, expected_output)

        # check the contents of the file
        output_file = expected_output['output_file']['path']
        with open(output_file) as fin:
            reader = csv.DictReader(fin)
            rows = [ row for row in reader ]

        self.assertEqual(len(rows), 2)
        self.assertDictEqual(rows[0], OrderedDict([('seg.mean', '3.1416')]))
        self.assertDictEqual(rows[1], OrderedDict([('seg.mean', '2.7183')]))


if __name__ == "__main__":
    unittest.main()
