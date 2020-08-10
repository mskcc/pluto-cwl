#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the run-facets-wrapper.cwl file
"""
import os
import json
import unittest
from tempfile import TemporaryDirectory, NamedTemporaryFile


# relative imports, from CLI and from parent project
if __name__ != "__main__":
    from .tools import run_command
    from .settings import CWL_DIR, CWL_ARGS, DATA_SETS

if __name__ == "__main__":
    from tools import run_command
    from settings import CWL_DIR, CWL_ARGS, DATA_SETS

cwl_file = os.path.join(CWL_DIR, 'run-facets-legacy-wrapper.cwl')

class TestRunFacetsWrapperCWL(unittest.TestCase):
    def test_run_facets_wrapper(self):
        """
        """
        input_json = {
            "snp_pileup": {
                "path": os.path.join(DATA_SETS['Proj_08390_G']['SNP_PILEUP_DIR'], "Sample24.Sample23.snp_pileup.gz"),
                "class": "File"
            },
            "sample_id": "Sample24"
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
                'arm_level_txt': {
                    'location': 'file://' + os.path.join(output_dir, 'Sample24.arm_level.txt'),
                    'basename': 'Sample24.arm_level.txt',
                    'class': 'File',
                    'checksum': 'sha1$df37c54ae4969257e436a7a7a595c42ef19ecbb5',
                    'size': 1824,
                    'path': os.path.join(output_dir, 'Sample24.arm_level.txt')
                    },
                'failed_txt': None,
                'gene_level_txt': {
                    'location': 'file://' + os.path.join(output_dir, 'Sample24.gene_level.txt'),
                    'basename': 'Sample24.gene_level.txt',
                    'class': 'File',
                    'checksum': 'sha1$4e916a52458151007486bf536acfff539fdc2ecc',
                    'size': 148195,
                    'path': os.path.join(output_dir, 'Sample24.gene_level.txt')
                },
                'hisens_Rdata': {
                    'location': 'file://' + os.path.join(output_dir, 'Sample24_hisens.Rdata'),
                    'basename': 'Sample24_hisens.Rdata',
                    'class': 'File',
                    'checksum': 'sha1$6cc3fd1fad17111e32c7c88b259a523092539181',
                    'size': 214260,
                    'path': os.path.join(output_dir, 'Sample24_hisens.Rdata')
                },
                'hisens_seg': {
                    'location': 'file://' + os.path.join(output_dir, 'Sample24_hisens.seg'),
                    'basename': 'Sample24_hisens.seg',
                    'class': 'File',
                    'checksum': 'sha1$652f9c6d0245af49bac6ca67a089af7d4e46801b',
                    'size': 1897,
                    'path': os.path.join(output_dir, 'Sample24_hisens.seg')
                },
                'purity_Rdata': {
                    'location': 'file://' + os.path.join(output_dir, 'Sample24_purity.Rdata'),
                    'basename': 'Sample24_purity.Rdata',
                    'class': 'File',
                    'checksum': 'sha1$ca2c46ceebbb960a02fb960134ba1711983e71c8',
                    'size': 213542,
                    'path': os.path.join(output_dir, 'Sample24_purity.Rdata')
                },
                'purity_seg': {
                    'location': 'file://' + os.path.join(output_dir, 'Sample24_purity.seg'),
                    'basename': 'Sample24_purity.seg',
                    'class': 'File',
                    'checksum': 'sha1$591e6d8b432e1e910fe4fb4b1814508131f960c9',
                    'size': 1285,
                    'path': os.path.join(output_dir, 'Sample24_purity.seg')
                },
                'qc_txt': {
                    'location': 'file://' + os.path.join(output_dir, 'Sample24.qc.txt'),
                    'basename': 'Sample24.qc.txt',
                    'class': 'File',
                    'checksum': 'sha1$d4a36726a5fcb7b268aae02d97ce4e382e42d9f6',
                    'size': 1339,
                    'path': os.path.join(output_dir, 'Sample24.qc.txt')
                },
                'hisens_cncf_txt': {
                    'basename': 'Sample24_hisens.cncf.txt',
                    'checksum': 'sha1$db9131a33889a1cac82e3bd6b3f0e5e182c65105',
                    'class': 'File',
                    'location': 'file://' + os.path.join(output_dir, 'Sample24_hisens.cncf.txt'),
                    'path': os.path.join(output_dir, 'Sample24_hisens.cncf.txt'),
                    'size': 5238
                },
                'purity_cncf_txt': {
                    'basename': 'Sample24_purity.cncf.txt',
                    'checksum': 'sha1$b331530e1e46b5ba1bdcedeb67f2aa82da6ebc5f',
                    'class': 'File',
                    'location': 'file://' + os.path.join(output_dir, 'Sample24_purity.cncf.txt'),
                    'path': os.path.join(output_dir, 'Sample24_purity.cncf.txt'),
                    'size': 3630
                },
                'stderr_txt': {
                    'basename': 'facets_legacy_stderr.txt',
                    'checksum': 'sha1$c85edeaaf165ce25f609fa425a24164c2439c784',
                    'class': 'File',
                    'location': 'file://' + os.path.join(output_dir, 'facets_legacy_stderr.txt'),
                    'path': os.path.join(output_dir, 'facets_legacy_stderr.txt'),
                    'size': 142
                },
                'stdout_txt': {
                    'basename': 'facets_legacy_stdout.txt',
                    'checksum': 'sha1$b9d9dcc1e039c79e1ee79a397ac8d37e8aa47fa6',
                    'class': 'File',
                    'location': 'file://' + os.path.join(output_dir, 'facets_legacy_stdout.txt'),
                    'path': os.path.join(output_dir, 'facets_legacy_stdout.txt'),
                    'size': 77}
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

            with open(os.path.join(output_dir, 'Sample24_hisens.cncf.txt')) as fin:
                self.assertEqual(len(fin.readlines()), 37)

if __name__ == "__main__":
    unittest.main()
