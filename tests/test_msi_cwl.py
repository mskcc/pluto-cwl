#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the snp-pileup-wrapper.cwl file
"""
import os
import json
import unittest
from tempfile import TemporaryDirectory, NamedTemporaryFile


# relative imports, from CLI and from parent project
if __name__ != "__main__":
    from .tools import run_command
    from .settings import CWL_DIR, CWL_ARGS, FACETS_SNPS_VCF, DATA_SETS

if __name__ == "__main__":
    from tools import run_command
    from settings import CWL_DIR, CWL_ARGS, FACETS_SNPS_VCF, DATA_SETS

cwl_file = os.path.join('./', 'msi.cwl') #CWL_DIR

class MSI_CWL(unittest.TestCase):
    def test_msi_cwl(self):
        """
        """
        input_json = {
              "d":    {"class": "File", "path": "/work/ci/resources/request_files/msisensor/b37_known_somatic_microsatellites.list" },
              "n":    {"class": "File", "path": os.path.join(DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample23.rg.md.abra.printreads.bam") },
              "t":    {"class": "File", "path": os.path.join(DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample24.rg.md.abra.printreads.bam") },
              "o": "Sample24.Sample23"
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
                'output': {
                    'location': 'file://' + os.path.join(output_dir, "Sample24.Sample23"),
                    'basename': "Sample24.Sample23",
                    'class': 'File',
                    'checksum': 'sha1$576fa40d02976be7527b394f482d16768874b174',
                    'size': 61,
                    'path': os.path.join(output_dir, "Sample24.Sample23")
                    }
                }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

if __name__ == "__main__":
    unittest.main()
