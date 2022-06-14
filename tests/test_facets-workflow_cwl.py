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

        self.maxDiff = None
        expected_output = {
        'annotated_maf': [
        OFile(name = "Sample24.Sample23_hisens.ccf.portal.maf", hash = "d91a8e15c66429b09f1b7db41bc38bdfa0b84c64", size = 11996607, dir = os.path.join(output_dir, 'facets/Sample24.Sample23'))
        ],
        'arm_level_txt': [
        OFile(name = "Sample24.arm_level.txt", hash = "df37c54ae4969257e436a7a7a595c42ef19ecbb5", size = 1824, dir = os.path.join(output_dir, 'facets/Sample24.Sample23'))
        ],
        'facets_txt': [
        OFile(name = "Sample24.txt", hash = "a0fb3df832efc18a66a8a54e5609666da5f4d7d7", size = 529, dir = os.path.join(output_dir, 'facets/Sample24.Sample23'))
        ],
        'failed_pairs': [],
        'gene_level_txt': [
        OFile(name = "Sample24.gene_level.txt", hash = "4e916a52458151007486bf536acfff539fdc2ecc", size = 148195, dir = os.path.join(output_dir, 'facets/Sample24.Sample23'))
        ],
        'hisens_rds': [
        OFile(name = "Sample24_hisens.rds", hash = "6bfd6c7f29c49ec8ef538dd468a3b4626b05bda2", size = 213986, dir = os.path.join(output_dir, 'facets/Sample24.Sample23'))
        ],
        'output_dir': {
            'basename': 'facets',
            'class': 'Directory',
            'listing': [
                {
                    'basename': 'Sample24.Sample23',
                    'class': 'Directory',
                    'listing': [
                        {
                            'basename': 'Sample24.Sample23_hisens.ccf.portal.maf',
                            'checksum': 'sha1$d91a8e15c66429b09f1b7db41bc38bdfa0b84c64',
                            'class': 'File',
                            'location': 'file://' + os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24.Sample23_hisens.ccf.portal.maf'),
                            'path': os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24.Sample23_hisens.ccf.portal.maf'),
                            'size': 11996607
                        },
                        {
                            'basename': 'Sample24.arm_level.txt',
                            'checksum': 'sha1$df37c54ae4969257e436a7a7a595c42ef19ecbb5',
                            'class': 'File',
                            'location': 'file://' + os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24.arm_level.txt'),
                            'path': os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24.arm_level.txt'),
                            'size': 1824
                        },
                        {
                            'basename': 'Sample24.txt',
                            'checksum': 'sha1$a0fb3df832efc18a66a8a54e5609666da5f4d7d7',
                            'class': 'File',
                            'location': 'file://' + os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24.txt'),
                            'path': os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24.txt'),
                            'size': 529
                        },
                        {
                            'basename': 'Sample24.gene_level.txt',
                            'checksum': 'sha1$4e916a52458151007486bf536acfff539fdc2ecc',
                            'class': 'File',
                            'location': 'file://' + os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24.gene_level.txt'),
                            'path': os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24.gene_level.txt'),
                            'size': 148195
                        },
                        {
                            'basename': 'Sample24_hisens.cncf.txt',
                            'checksum': 'sha1$db9131a33889a1cac82e3bd6b3f0e5e182c65105',
                            'class': 'File',
                            'location': 'file://' + os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24_hisens.cncf.txt'),
                            'path': os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24_hisens.cncf.txt'),
                            'size': 5238
                        },
                        {
                            'basename': 'Sample24_hisens.rds',
                            'checksum': 'sha1$6bfd6c7f29c49ec8ef538dd468a3b4626b05bda2',
                            'class': 'File',
                            'location': 'file://' + os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24_hisens.rds'),
                            'path': os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24_hisens.rds'),
                            'size': 213986
                        },
                        {
                            'basename': 'Sample24_hisens.seg',
                            'checksum': 'sha1$652f9c6d0245af49bac6ca67a089af7d4e46801b',
                            'class': 'File',
                            'location': 'file://' + os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24_hisens.seg'),
                            'path': os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24_hisens.seg'),
                            'size': 1897
                        },
                        {
                            'basename': 'Sample24_purity.rds',
                            'checksum': 'sha1$dd8b967f84b191ff76214c6110db8d0e65f6514c',
                            'class': 'File',
                            'location': 'file://' + os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24_purity.rds'),
                            'path': os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24_purity.rds'),
                            'size': 213356
                        },
                        {
                            'basename': 'Sample24_purity.seg',
                            'checksum': 'sha1$591e6d8b432e1e910fe4fb4b1814508131f960c9',
                            'class': 'File',
                            'location': 'file://' + os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24_purity.seg'),
                            'path': os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24_purity.seg'),
                            'size': 1285
                        },
                        {
                            'basename': 'Sample24.qc.txt',
                            'checksum': 'sha1$d4a36726a5fcb7b268aae02d97ce4e382e42d9f6',
                            'class': 'File',
                            'location': 'file://' + os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24.qc.txt'),
                            'path': os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24.qc.txt'),
                            'size': 1339}
                    ],
                    'location': 'file://' + os.path.join(output_dir, 'facets/Sample24.Sample23'),
                    'path': os.path.join(output_dir, 'facets/Sample24.Sample23')
                }
            ],
            'location': 'file://' + os.path.join(output_dir, 'facets'),
            'path': os.path.join(output_dir, 'facets')
        },
        'hisens_seg': [{
            'location': 'file://' + os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24_hisens.seg'),
            'basename': 'Sample24_hisens.seg',
            'class': 'File',
            'checksum': 'sha1$652f9c6d0245af49bac6ca67a089af7d4e46801b',
            'size': 1897,
            'path': os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24_hisens.seg')
        }],
        'purity_rds': [{
            'location': 'file://' + os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24_purity.rds'),
            'basename': 'Sample24_purity.rds',
            'class': 'File',
            'checksum': 'sha1$dd8b967f84b191ff76214c6110db8d0e65f6514c',
            'size': 213356,
            'path': os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24_purity.rds')
        }],
        'purity_seg': [{
            'location': 'file://' + os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24_purity.seg'),
            'basename': 'Sample24_purity.seg',
            'class': 'File',
            'checksum': 'sha1$591e6d8b432e1e910fe4fb4b1814508131f960c9',
            'size': 1285,
            'path': os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24_purity.seg')
        }],
        'qc_txt': [{
            'location': 'file://' + os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24.qc.txt'),
            'basename': 'Sample24.qc.txt',
            'class': 'File',
            'checksum': 'sha1$d4a36726a5fcb7b268aae02d97ce4e382e42d9f6',
            'size': 1339,
            'path': os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24.qc.txt')
        }],
        'hisens_cncf_txt': [{
            'basename': 'Sample24_hisens.cncf.txt',
            'checksum': 'sha1$db9131a33889a1cac82e3bd6b3f0e5e182c65105',
            'class': 'File',
            'location': 'file://' + os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24_hisens.cncf.txt'),
            'path': os.path.join(output_dir, 'facets/Sample24.Sample23/Sample24_hisens.cncf.txt'),
            'size': 5238
           }],
        }

        self.assertCWLDictEqual(output_json, expected_output)

if __name__ == "__main__":
    unittest.main()
