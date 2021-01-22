#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unit tests for the fusion_filter.cwl
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
from pluto.settings import CWL_ARGS, DATA_SETS, KNOWN_FUSIONS_FILE
sys.path.pop(0)

cwl_file = CWLFile('fusion_filter.cwl')

class TestFusionFilter(unittest.TestCase):
    def test_fusion_filter1(self):
        """
        """
        fusion_file = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.portal.txt")

        input_json = {
            "fusions_file": {
                  "class": "File",
                  "path": fusion_file
                },
            "output_filename": "data_fusions.txt",
            "known_fusions_file": {
                "class": "File",
                "path": KNOWN_FUSIONS_FILE
            }
        }

        with TemporaryDirectory() as tmpdir:

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

            expected_output = {
                "output_file": {
                    'location': 'file://' + os.path.join(output_dir, "data_fusions.txt"),
                    'basename': "data_fusions.txt",
                    'class': 'File',
                    'checksum': 'sha1$c16f763b248813fcdde76f7486f1ddc4e9856038',
                    'size': 99,
                    'path': os.path.join(output_dir, "data_fusions.txt")
                    }
                }
        self.assertDictEqual(output_json, expected_output)




if __name__ == "__main__":
    unittest.main()
