#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the annotate-maf-wrapper.cwl file
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
from pluto.settings import CWL_ARGS, DATA_SETS
sys.path.pop(0)

cwl_file = CWLFile('annotate-maf-wrapper.cwl')

class TestAnnotateMafWrapperCWL(unittest.TestCase):
    def test_run_facets_wrapper(self):
        """
        """
        input_maf = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf")
        input_rds = os.path.join(DATA_SETS['Proj_08390_G']['FACETS_SUITE_DIR'], "Sample1_hisens.rds")
        input_json = {
            "maf_file": {
                "path": input_maf,
                "class": "File"
            },
            "facets_rds": {
                "path": input_rds,
                "class": "File"
            },
            "output_filename": "Sample1_hisens.ccf.maf"
        }
        with TemporaryDirectory() as tmpdir:
            input_json_file = os.path.join(tmpdir, "input.json")
            with open(input_json_file, "w") as json_out:
                json.dump(input_json, json_out)

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
                'failed_txt': None,
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'Sample1_hisens.ccf.maf'),
                    'basename': 'Sample1_hisens.ccf.maf',
                    'class': 'File',
                    'checksum': 'sha1$7e478a8a44d27735f26e368989c672ed6ef5d52a',
                    'size': 19217199,
                    'path': os.path.join(output_dir, 'Sample1_hisens.ccf.maf')
                },
                'stderr_txt': {
                    'basename': 'annotate_maf_stderr.txt',
                    'checksum': 'sha1$2e672f99c23a2d827c1d33e06377870cdd9c8090',
                    'class': 'File',
                    'location': 'file://' + os.path.join(output_dir,'annotate_maf_stderr.txt'),
                    'path': os.path.join(output_dir,'annotate_maf_stderr.txt'),
                    'size': 105
                },
               'stdout_txt': {
                    'basename': 'annotate_maf_stdout.txt',
                    'checksum': 'sha1$da39a3ee5e6b4b0d3255bfef95601890afd80709',
                    'class': 'File',
                    'location': 'file://' + os.path.join(output_dir,'annotate_maf_stdout.txt'),
                    'path': os.path.join(output_dir,'annotate_maf_stdout.txt'),
                    'size': 0
                }
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

            # with open(output_json['stdout_txt']["path"]) as f:
            #     for line in f:
            #         print(line)
            #
            # with open(output_json['stderr_txt']["path"]) as f:
            #     for line in f:
            #         print(line)

if __name__ == "__main__":
    unittest.main()
