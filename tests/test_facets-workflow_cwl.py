#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unit tests for the facets-workflow.cwl
"""
import os
import json
import unittest
from tempfile import TemporaryDirectory, NamedTemporaryFile


# relative imports, from CLI and from parent project
if __name__ != "__main__":
    from .tools import run_command
    from .settings import CWL_DIR, CWL_ARGS, DATA_SETS, FACETS_SNPS_VCF

if __name__ == "__main__":
    from tools import run_command
    from settings import CWL_DIR, CWL_ARGS, DATA_SETS, FACETS_SNPS_VCF

cwl_file = os.path.join(CWL_DIR, 'facets-workflow.cwl')

class TestFacetsWorkflow(unittest.TestCase):
    def test_facets_workflow(self):
        """
        """
        # the smallest pair of bam files in the test dataset
        tumor_bam = os.path.join(DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample24.rg.md.abra.printreads.bam")
        normal_bam = os.path.join(DATA_SETS['Proj_08390_G']['BAM_DIR'], "Sample23.rg.md.abra.printreads.bam")
        pair_maf = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample24.Sample23.muts.maf")
        input_json = {
            "snps_vcf": {
                "path": FACETS_SNPS_VCF,
                "class": "File"
            },
            "pairs": [
                {
                    "normal_bam": {
                        "path": normal_bam,
                        "class": "File"
                    },
                    "tumor_bam": {
                        "path": tumor_bam,
                        "class": "File"
                    },
                    "pair_maf": {
                        "path": pair_maf,
                        "class": "File"
                    },
                    "pair_id": "Sample24.Sample23",
                    "normal_id": "Sample23",
                    "tumor_id": "Sample24"
                }
            ]
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

            expected_output = {'annotated_maf': [{'basename': 'Sample24.Sample23_hisens.ccf.portal.maf',
                     'checksum': 'sha1$d91a8e15c66429b09f1b7db41bc38bdfa0b84c64',
                     'class': 'File',
                     'location': 'file://' + os.path.join(output_dir, 'Sample24.Sample23_hisens.ccf.portal.maf'),
                     'path': os.path.join(output_dir, 'Sample24.Sample23_hisens.ccf.portal.maf'),
                     'size': 11996607}],
  'arm_level_txt': [{'basename': 'Sample24.Sample23.arm_level.txt',
                     'checksum': 'sha1$9f89c90a0f196ab2771f5e9aadc5893a5e509609',
                     'class': 'File',
                     'location': 'file://' + os.path.join(output_dir, 'Sample24.Sample23.arm_level.txt'),
                     'path': os.path.join(output_dir, 'Sample24.Sample23.arm_level.txt'),
                     'size': 2148}],
  'facets_txt': [{'basename': 'Sample24.Sample23.txt',
                  'checksum': 'sha1$b1cc75293a4bddfe36540a6ea4248a44acae6db0',
                  'class': 'File',
                  'location': 'file://' + os.path.join(output_dir, 'Sample24.Sample23.txt'),
                  'path': os.path.join(output_dir, 'Sample24.Sample23.txt'),
                  'size': 498}],
  'gene_level_txt': [{'basename': 'Sample24.Sample23.gene_level.txt',
                      'checksum': 'sha1$f056daac29597163fecb90e69d801953a74c7619',
                      'class': 'File',
                      'location': 'file://' + os.path.join(output_dir, 'Sample24.Sample23.gene_level.txt'),
                      'path': os.path.join(output_dir, 'Sample24.Sample23.gene_level.txt'),
                      'size': 156007}],
  'hisens_rds': [{'basename': 'Sample24.Sample23_hisens.rds',
                  'checksum': 'sha1$6bfd6c7f29c49ec8ef538dd468a3b4626b05bda2',
                  'class': 'File',
                  'location': 'file://' + os.path.join(output_dir, 'Sample24.Sample23_hisens.rds'),
                  'path': os.path.join(output_dir, 'Sample24.Sample23_hisens.rds'),
                  'size': 213986}],
  'hisens_seg': [{'basename': 'Sample24.Sample23_hisens.seg',
                  'checksum': 'sha1$1ddb1a4a43c5707188ae597792c372eee2910d1c',
                  'class': 'File',
                  'location': 'file://' + os.path.join(output_dir, 'Sample24.Sample23_hisens.seg'),
                  'path': os.path.join(output_dir, 'Sample24.Sample23_hisens.seg'),
                  'size': 2221}],
  'purity_rds': [{'basename': 'Sample24.Sample23_purity.rds',
                  'checksum': 'sha1$dd8b967f84b191ff76214c6110db8d0e65f6514c',
                  'class': 'File',
                  'location': 'file://' + os.path.join(output_dir, 'Sample24.Sample23_purity.rds'),
                  'path': os.path.join(output_dir, 'Sample24.Sample23_purity.rds'),
                  'size': 213356}],
  'purity_seg': [{'basename': 'Sample24.Sample23_purity.seg',
                  'checksum': 'sha1$dadb691707754772520a82c9148620d428281328',
                  'class': 'File',
                  'location': 'file://' + os.path.join(output_dir, 'Sample24.Sample23_purity.seg'),
                  'path': os.path.join(output_dir, 'Sample24.Sample23_purity.seg'),
                  'size': 1501}],
  'qc_txt': [{'basename': 'Sample24.Sample23.qc.txt',
              'checksum': 'sha1$704180fe9d6ccf38a7332d40b524bea3a2d3ff98',
              'class': 'File',
              'location': 'file://' + os.path.join(output_dir, 'Sample24.Sample23.qc.txt'),
              'path': os.path.join(output_dir, 'Sample24.Sample23.qc.txt'),
              'size': 1357}],
  'snp_pileup': [{'basename': 'Sample24.Sample23.snp_pileup.gz',
                  'checksum': 'sha1$755a8b64f45c819b4e2c481e64bf2fe36d1f5361',
                  'class': 'File',
                  'location': 'file://' + os.path.join(output_dir, 'Sample24.Sample23.snp_pileup.gz'),
                  'path': os.path.join(output_dir, 'Sample24.Sample23.snp_pileup.gz'),
                  'size': 34851004}]}

            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

if __name__ == "__main__":
    unittest.main()
