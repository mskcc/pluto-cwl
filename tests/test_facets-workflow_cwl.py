#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unit tests for the facets-workflow.cwl
"""
import os
import sys
from datasets import (
    DATA_SETS,
)
from pluto import (
    PlutoTestCase,
    CWLFile,
    OFile,
    ODir
)


class TestFacetsWorkflow(PlutoTestCase):
    cwl_file = CWLFile('facets-workflow.cwl')

    def test_facets_workflow(self):
        """
        """
        # the smallest pair of bam files in the test dataset
        pair_maf24 = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample24.Sample23.muts.maf")
        snp_pileup24 = os.path.join(DATA_SETS['Proj_08390_G']['SNP_PILEUP_DIR'], "Sample24.Sample23.snp_pileup.gz")

        pair_maf12 = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample12.Sample11.muts.maf")
        snp_pileup12 = os.path.join(DATA_SETS['Proj_08390_G']['SNP_PILEUP_DIR'], "Sample12.Sample11.snp_pileup.gz")

        self.input = {
            "pairs": [
                {
                    "pair_maf": {"path": pair_maf24, "class": "File"},
                    "snp_pileup": {"path": snp_pileup24, "class": "File"},
                    "pair_id": "Sample24.Sample23",
                    "normal_id": "Sample23",
                    "tumor_id": "Sample24"
                },
                {
                    "pair_maf": {"path": pair_maf12, "class": "File"},
                    "snp_pileup": {"path": snp_pileup12, "class": "File"},
                    "pair_id": "Sample12.Sample11",
                    "normal_id": "Sample11",
                    "tumor_id": "Sample12"
                },
            ]
        }

        output_json, output_dir = self.run_cwl()
        # output_dir_sample24 = os.path.join(output_dir, 'facets/Sample24.Sample23')
        # output_dir_sample24 = os.path.join(output_dir, 'Sample24.Sample23')

        self.maxDiff = None
        expected_output = {
        "pairs":[{
            "pair_id": "Sample24.Sample23",
            "normal_id": "Sample23",
            "tumor_id": "Sample24",
            'annotated_maf': OFile(name = "Sample24.Sample23_hisens.ccf.portal.maf", hash = "d91a8e15c66429b09f1b7db41bc38bdfa0b84c64", size = 11996607, dir = output_dir),
            'arm_level_txt': OFile(name = "Sample24.arm_level.txt", hash = "df37c54ae4969257e436a7a7a595c42ef19ecbb5", size = 1824, dir = output_dir),
            'facets_txt': OFile(name = "Sample24.txt", hash = "a0fb3df832efc18a66a8a54e5609666da5f4d7d7", size = 529, dir = output_dir),
            'gene_level_txt': OFile(name = "Sample24.gene_level.txt", hash = "4e916a52458151007486bf536acfff539fdc2ecc", size = 148195, dir = output_dir),
            'hisens_rds': OFile(name = "Sample24_hisens.rds", hash = "6bfd6c7f29c49ec8ef538dd468a3b4626b05bda2", size = 213986, dir = output_dir),
            'hisens_png': OFile(name = 'Sample24_hisens.png', size = 168166, hash = '6af56798d0d8e3b49c26ab9d0adc855c3c8a5a50', dir = output_dir),
            'hisens_seg': OFile(name = "Sample24_hisens.seg", hash = "652f9c6d0245af49bac6ca67a089af7d4e46801b", size = 1897, dir = output_dir),
            'purity_png': OFile(name = 'Sample24_purity.png', size = 164021, hash = '7db765d900c8a431ab0325098b81eda2cd0780bf', dir = output_dir),
            'purity_rds': OFile(name = "Sample24_purity.rds", hash = "dd8b967f84b191ff76214c6110db8d0e65f6514c", size = 213356, dir = output_dir),
            'purity_seg': OFile(name = "Sample24_purity.seg", hash = "591e6d8b432e1e910fe4fb4b1814508131f960c9", size = 1285, dir = output_dir),
            'qc_txt': OFile(name = "Sample24.qc.txt", hash = "d4a36726a5fcb7b268aae02d97ce4e382e42d9f6", size = 1339, dir = output_dir),
            'hisens_cncf_txt': OFile(name = "Sample24_hisens.cncf.txt", hash = "db9131a33889a1cac82e3bd6b3f0e5e182c65105", size = 5238, dir = output_dir)
        },
        {
            "pair_id": "Sample12.Sample11",
            "normal_id": "Sample11",
            "tumor_id": "Sample12",
            'annotated_maf': OFile(name = "Sample12.Sample11_hisens.ccf.portal.maf", hash = "6b570289519f3c5d2b2d52cd3f4ce7954f8d083f", size = 19857901, dir = output_dir),
            'arm_level_txt': OFile(name = "Sample12.arm_level.txt", hash = "d53fc0f25763c9171a254af536f28dbe476a093f", size = 1453, dir = output_dir),
            'facets_txt': OFile(name = "Sample12.txt", hash = "bc2449d989b00bf3dc7ac4dff3b8fcfeeea2f84e", size = 467, dir = output_dir),
            'gene_level_txt': OFile(name = "Sample12.gene_level.txt", hash = "9f5b115b1499b59f319b1a69035de32def09242c", size = 176045, dir = output_dir),
            'hisens_rds': OFile(name = "Sample12_hisens.rds", hash = "8d20a8e4393ba858c65962f870924d476d73763a", size = 303751, dir = output_dir),
            'hisens_png': OFile(name = 'Sample12_hisens.png', size = 163006, hash = '31d748a89932682291e1f4c05f21a98e50ed38bf', dir = output_dir),
            'hisens_seg': OFile(name = "Sample12_hisens.seg", hash = "efe97bf22666b58ea61e026d681b3007af40319f", size = 2907, dir = output_dir),
            'purity_png': OFile(name = 'Sample12_purity.png', size = 158950, hash = '77cbac10e777fa8f2d0f981e04af0eab8efb00cc', dir = output_dir),
            'purity_rds': OFile(name = "Sample12_purity.rds", hash = "9c2bd66081e4cd458f4407c98a748a7336388212", size = 302766, dir = output_dir),
            'purity_seg': OFile(name = "Sample12_purity.seg", hash = "4c47f0c5b7e002aa9e745b4e31d68c948fa7710d", size = 2079, dir = output_dir),
            'qc_txt': OFile(name = "Sample12.qc.txt", hash = "827ec19a3d567415c9f69237e796cda1f7323c90", size = 1323, dir = output_dir),
            'hisens_cncf_txt': OFile(name = "Sample12_hisens.cncf.txt", hash = "ef857c6242c12b4d4e0cccd7f1eb2df7fd2bbc7b", size = 8207, dir = output_dir)
        }]
        }

        self.assertCWLDictEqual(output_json, expected_output)
