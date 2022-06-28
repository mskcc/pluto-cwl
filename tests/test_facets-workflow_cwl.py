#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unit tests for the facets-workflow.cwl
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile
from pluto.settings import DATA_SETS
from pluto.serializer import OFile, ODir
sys.path.pop(0)

class TestFacetsWorkflow(PlutoTestCase):
    cwl_file = CWLFile('facets-workflow.cwl')

    def test_facets_workflow(self):
        """
        """
        # the smallest pair of bam files in the test dataset
        pair_maf = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample24.Sample23.muts.maf")
        snp_pileup = os.path.join(DATA_SETS['Proj_08390_G']['SNP_PILEUP_DIR'], "Sample24.Sample23.snp_pileup.gz")
        self.input = {
            "pairs": [
                {
                    "pair_maf": {
                        "path": pair_maf,
                        "class": "File"
                    },
                    "snp_pileup": {
                        "path": snp_pileup,
                        "class": "File"
                    },
                    "pair_id": "Sample24.Sample23",
                    "normal_id": "Sample23",
                    "tumor_id": "Sample24"
                }
            ]
        }

        output_json, output_dir = self.run_cwl()
        output_dir_sample24 = os.path.join(output_dir, 'facets/Sample24.Sample23')

        self.maxDiff = None
        expected_output = {
        'annotated_maf': [
            OFile(name = "Sample24.Sample23_hisens.ccf.portal.maf", hash = "d91a8e15c66429b09f1b7db41bc38bdfa0b84c64", size = 11996607, dir = output_dir_sample24)
            ],
        'arm_level_txt': [
            OFile(name = "Sample24.arm_level.txt", hash = "df37c54ae4969257e436a7a7a595c42ef19ecbb5", size = 1824, dir = output_dir_sample24)
            ],
        'facets_txt': [
            OFile(name = "Sample24.txt", hash = "a0fb3df832efc18a66a8a54e5609666da5f4d7d7", size = 529, dir = output_dir_sample24)
            ],
        'failed_pairs': [],
        'gene_level_txt': [
            OFile(name = "Sample24.gene_level.txt", hash = "4e916a52458151007486bf536acfff539fdc2ecc", size = 148195, dir = output_dir_sample24)
            ],
        'hisens_rds': [
            OFile(name = "Sample24_hisens.rds", hash = "6bfd6c7f29c49ec8ef538dd468a3b4626b05bda2", size = 213986, dir = output_dir_sample24)
            ],
        'output_dir': ODir(name = 'facets', dir = output_dir, items = [
            ODir(name = 'Sample24.Sample23', items = [
                OFile(name = 'Sample24.Sample23_hisens.ccf.portal.maf', size = 11996607, hash = 'd91a8e15c66429b09f1b7db41bc38bdfa0b84c64'),
                OFile(name = 'Sample24.arm_level.txt', size = 1824, hash = 'df37c54ae4969257e436a7a7a595c42ef19ecbb5'),
                OFile(name = 'Sample24.txt', size = 529, hash = 'a0fb3df832efc18a66a8a54e5609666da5f4d7d7'),
                OFile(name = 'Sample24.gene_level.txt', size = 148195, hash = '4e916a52458151007486bf536acfff539fdc2ecc'),
                OFile(name = 'Sample24_hisens.cncf.txt', size = 5238, hash = 'db9131a33889a1cac82e3bd6b3f0e5e182c65105'),
                OFile(name = 'Sample24_hisens.rds', size = 213986, hash = '6bfd6c7f29c49ec8ef538dd468a3b4626b05bda2'),
                OFile(name = 'Sample24_hisens.seg', size = 1897, hash = '652f9c6d0245af49bac6ca67a089af7d4e46801b'),
                OFile(name = 'Sample24_purity.rds', size = 213356, hash = 'dd8b967f84b191ff76214c6110db8d0e65f6514c'),
                OFile(name = 'Sample24_purity.seg', size = 1285, hash = '591e6d8b432e1e910fe4fb4b1814508131f960c9'),
                OFile(name = 'Sample24.qc.txt', size = 1339, hash = 'd4a36726a5fcb7b268aae02d97ce4e382e42d9f6'),
                ])
            ]),
        'hisens_seg': [
            OFile(name = "Sample24_hisens.seg", hash = "652f9c6d0245af49bac6ca67a089af7d4e46801b", size = 1897, dir = output_dir_sample24)
            ],
        'purity_rds': [
            OFile(name = "Sample24_purity.rds", hash = "dd8b967f84b191ff76214c6110db8d0e65f6514c", size = 213356, dir = output_dir_sample24)
            ],
        'purity_seg': [
            OFile(name = "Sample24_purity.seg", hash = "591e6d8b432e1e910fe4fb4b1814508131f960c9", size = 1285, dir = output_dir_sample24)
            ],
        'qc_txt': [
            OFile(name = "Sample24.qc.txt", hash = "d4a36726a5fcb7b268aae02d97ce4e382e42d9f6", size = 1339, dir = output_dir_sample24)
            ],
        'hisens_cncf_txt': [
            OFile(name = "Sample24_hisens.cncf.txt", hash = "db9131a33889a1cac82e3bd6b3f0e5e182c65105", size = 5238, dir = output_dir_sample24)
               ],
        }

        self.assertCWLDictEqual(output_json, expected_output)

if __name__ == "__main__":
    unittest.main()
