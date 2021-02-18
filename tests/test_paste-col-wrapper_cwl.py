#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the paste-col-wrapper.cwl
"""
import os
import sys
import json
import unittest
from tempfile import TemporaryDirectory

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import run_command, CWLFile
from pluto.settings import CWL_ARGS
sys.path.pop(0)

cwl_file = CWLFile('paste-col-wrapper.cwl')

class TestPasteCol(unittest.TestCase):
    def test_paste_col_1(self):
        """
        """
        with TemporaryDirectory() as tmpdir:
            # make a dummy file with some lines
            input_lines = ["HEADER1", "foo1", "bar1"]
            input_file = os.path.join(tmpdir, "input.txt")
            with open(input_file, "w") as fout:
                for line in input_lines:
                    fout.write(line + '\n')

            input_json = {
                "input_file": {
                      "class": "File",
                      "path": input_file
                    },
                "output_filename": "output.txt",
                "header": "HEADER2",
                "value": "foo2"
                }
            input_json_file = os.path.join(tmpdir, "input.json")
            with open(input_json_file, "w") as input_json_file_data:
                json.dump(input_json, input_json_file_data)

            output_dir = os.path.join(tmpdir, "output")
            tmp_dir = os.path.join(tmpdir, "tmp")
            cache_dir = os.path.join(tmpdir, "cache")

            command = [
            "cwl-runner",
            *CWL_ARGS,
            "--outdir", output_dir,
            "--tmpdir-prefix", tmp_dir,
            "--cachedir", cache_dir,
            cwl_file, input_json_file
            ]

            returncode, proc_stdout, proc_stderr = run_command(command)

            if returncode != 0:
                print(proc_stderr)

            self.assertEqual(returncode, 0)

            output_json = json.loads(proc_stdout)

            # check the contents of the concatenated file; should be the same as the input
            output_file = output_json['output_file']['path']
            with open(output_file) as fin:
                output_lines = [ line.strip() for line in fin ]

            expected_lines = ['HEADER1\tHEADER2', 'foo1\tfoo2', 'bar1\tfoo2']
            self.assertEqual(output_lines, expected_lines)

            expected_output = {
                'failed_txt': None,
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'output.txt'),
                    'basename': 'output.txt',
                    'class': 'File',
                    'checksum': 'sha1$34753fd98b2355d54740f3fdfc6490262c15dd59',
                    'size': 36,
                    'path': os.path.join(output_dir, 'output.txt')
                    },
                'stderr_txt': {
                    'basename': 'output.txt_stderr.txt',
                    'checksum': 'sha1$da39a3ee5e6b4b0d3255bfef95601890afd80709',
                    'class': 'File',
                    'location': 'file://' + os.path.join(output_dir, 'output.txt_stderr.txt'),
                    'path': os.path.join(output_dir, 'output.txt_stderr.txt'),
                    'size': 0
                },
                'stdout_txt': {
                    'basename': 'output.txt_stdout.txt',
                    'checksum': 'sha1$da39a3ee5e6b4b0d3255bfef95601890afd80709',
                    'class': 'File',
                    'location': 'file://' + os.path.join(output_dir, 'output.txt_stdout.txt'),
                    'path': os.path.join(output_dir, 'output.txt_stdout.txt'),
                    'size': 0
                }
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)



if __name__ == "__main__":
    unittest.main()
