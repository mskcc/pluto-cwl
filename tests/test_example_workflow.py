#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unit tests for the example_workflow.cwl
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

class TestExampleWorkflow(PlutoTestCase):
    cwl_file = CWLFile('example_workflow.cwl')

    def test_example_workflow(self):
        """
        Test case for the example workflow
        """
        self.maxDiff = None
        self.input = {
            'value': "ABC",
            "samples": [
                {"sample_id": "1"},
                {"sample_id": "2"}
            ]
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            "output_file": OFile(name = 'output.concat.tsv', hash = 'd4297dfdad25ac92ffae2ce61c6cfe12c4089c28', size = 27, dir = output_dir),
            'env': OFile(name = 'env.txt', dir = output_dir)
        }
        strip_related_keys = [
        ('basename', 'env.txt', ['size', 'checksum']),
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)

        output_file = os.path.join(output_dir, "output.concat.tsv")
        with open(output_file) as f:
            lines = [ l.strip() for l in f ]

        expected_lines = [
        'SampleID\tValue',
        '1\tABC',
        '2\tABC',
        ]
        self.assertEqual(lines, expected_lines)


if __name__ == "__main__":
    unittest.main()
