#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the run-facets-wrapper.cwl file
"""
import os
import sys
import unittest



from pluto import (
    PlutoTestCase, 
    CWLFile,
    DATA_SETS,
    OFile
)


class TestRunFacetsWrapperCWL(PlutoTestCase):
    cwl_file = CWLFile('run-facets-wrapper.cwl')

    def test_run_facets_wrapper(self):
        """
        """
        self.input = {
            "snp_pileup": {
                "path": os.path.join(DATA_SETS['Proj_08390_G']['SNP_PILEUP_DIR'], "Sample24.Sample23.snp_pileup.gz"),
                "class": "File"
            },
            "sample_id": "Sample24"
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
        'failed': False,
        'arm_level_txt': OFile(name = "Sample24.arm_level.txt", hash = "df37c54ae4969257e436a7a7a595c42ef19ecbb5", size = 1824, dir = output_dir),
        'gene_level_txt': OFile(name = "Sample24.gene_level.txt", hash = "4e916a52458151007486bf536acfff539fdc2ecc", size = 148195, dir = output_dir),
        'hisens_rds': OFile(name = "Sample24_hisens.rds", hash = "6bfd6c7f29c49ec8ef538dd468a3b4626b05bda2", size = 213986, dir = output_dir),
        'hisens_seg': OFile(name = "Sample24_hisens.seg", hash = "652f9c6d0245af49bac6ca67a089af7d4e46801b", size = 1897, dir = output_dir),
        'hisens_png': OFile(name = "Sample24_hisens.png", hash = "6af56798d0d8e3b49c26ab9d0adc855c3c8a5a50", size = 168166, dir = output_dir),
        'purity_rds': OFile(name = "Sample24_purity.rds", hash = "dd8b967f84b191ff76214c6110db8d0e65f6514c", size = 213356, dir = output_dir),
        'purity_png': OFile(name = "Sample24_purity.png", hash = "7db765d900c8a431ab0325098b81eda2cd0780bf", size = 164021, dir = output_dir),
        'purity_seg': OFile(name = "Sample24_purity.seg", hash = "591e6d8b432e1e910fe4fb4b1814508131f960c9", size = 1285, dir = output_dir),
        'qc_txt': OFile(name = "Sample24.qc.txt", hash = "d4a36726a5fcb7b268aae02d97ce4e382e42d9f6", size = 1339, dir = output_dir),
        'output_txt': OFile(name = "Sample24.txt", hash = "4769dc7b8d4b127383e1936c07cdba1e2e09aecb", size = 480, dir = output_dir),
        "stderr_txt": OFile(name = "facets_stderr.txt", dir = output_dir),
        "stdout_txt": OFile(name = "facets_stdout.txt", dir = output_dir),
        }
        strip_related_keys = [
            ('basename', 'facets_stderr.txt', ['size', 'checksum']),
            ('basename', 'facets_stdout.txt', ['size', 'checksum']),
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)

        with open(os.path.join(output_dir, 'Sample24_hisens.seg')) as fin:
            self.assertEqual(len(fin.readlines()), 37)
        with open(os.path.join(output_dir, 'Sample24_purity.seg')) as fin:
            self.assertEqual(len(fin.readlines()), 25)



