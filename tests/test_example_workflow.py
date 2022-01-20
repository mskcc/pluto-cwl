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
            "output_file": {
                "location": "file://" + os.path.join(output_dir, "output.concat.tsv"),
                "basename": "output.concat.tsv",
                "class": "File",
                "checksum": "sha1$be6fb2e96f81c63a0b5fc6392a317ba3afbbca19",
                "size": 30,
                "path": os.path.join(output_dir, "output.concat.tsv")
            },
            'env': {
                'basename': 'env.txt',
                # 'checksum': 'sha1$e2f2bf6581461560dc9d4c4c970b5b7b1ba15852',
                'class': 'File',
                'location':  "file://" + os.path.join(output_dir, "env.txt"),
                'path': os.path.join(output_dir, 'env.txt')
                # 'size': 456
            }
        }
        output_json['env'].pop('checksum')
        output_json['env'].pop('size')

        self.assertCWLDictEqual(output_json, expected_output)

        output_file = output_json['output_file']['path']
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
