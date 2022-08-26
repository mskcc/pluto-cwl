#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the TMB analysis workflow cwl which uses multiple input samples
"""
import os
import sys
import unittest

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import PlutoTestCase, CWLFile
from pluto.serializer import OFile
sys.path.pop(0)


class TestTmbWorkflow(PlutoTestCase):
    cwl_file = CWLFile('tmb_workflow.cwl')

    def test_tmb_workflow(self):
        """
        Test case with multiple full size samples and IMPACT505 assay coverage value
        """
        self.input = {
            "assay_coverage":  '1213770', # IMPACT505
            "pairs": [
                {
                    "pair_id": "Sample10.Sample9", "tumor_id": "Sample10", "normal_id": "Sample9",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample10.Sample9.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample12.Sample11", "tumor_id": "Sample12", "normal_id": "Sample11",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample12.Sample11.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample14.Sample13", "tumor_id": "Sample14", "normal_id": "Sample13",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample14.Sample13.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample16.Sample15", "tumor_id": "Sample16", "normal_id": "Sample15",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample16.Sample15.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample18.Sample17", "tumor_id": "Sample18", "normal_id": "Sample17",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample18.Sample17.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample20.Sample19", "tumor_id": "Sample20", "normal_id": "Sample19",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample20.Sample19.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample22.Sample21", "tumor_id": "Sample22", "normal_id": "Sample21",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample22.Sample21.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample24.Sample23", "tumor_id": "Sample24", "normal_id": "Sample23",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample24.Sample23.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample26.Sample25", "tumor_id": "Sample26", "normal_id": "Sample25",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample26.Sample25.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample28.Sample27", "tumor_id": "Sample28", "normal_id": "Sample27",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample28.Sample27.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample30.Sample29", "tumor_id": "Sample30", "normal_id": "Sample29",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample30.Sample29.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample32.Sample31", "tumor_id": "Sample32", "normal_id": "Sample31",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample32.Sample31.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample34.Sample33", "tumor_id": "Sample34", "normal_id": "Sample33",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample34.Sample33.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample50.Sample49", "tumor_id": "Sample50", "normal_id": "Sample49",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample50.Sample49.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample52.Sample51", "tumor_id": "Sample52", "normal_id": "Sample51",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample52.Sample51.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample54.Sample53", "tumor_id": "Sample54", "normal_id": "Sample53",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample54.Sample53.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample62.Sample61", "tumor_id": "Sample62", "normal_id": "Sample61",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample62.Sample61.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample70.Sample69", "tumor_id": "Sample70", "normal_id": "Sample69",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample70.Sample69.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample74.Sample73", "tumor_id": "Sample74", "normal_id": "Sample73",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample74.Sample73.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample80.Sample79", "tumor_id": "Sample80", "normal_id": "Sample79",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample80.Sample79.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample84.Sample83", "tumor_id": "Sample84", "normal_id": "Sample83",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample84.Sample83.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample90.Sample89", "tumor_id": "Sample90", "normal_id": "Sample89",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample90.Sample89.muts.maf"), "class": "File"}
                },
                {
                    "pair_id": "Sample92.Sample91", "tumor_id": "Sample92", "normal_id": "Sample91",
                    "pair_maf": {"path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample92.Sample91.muts.maf"), "class": "File"}
                },
                ]
            }
        output_json, output_dir = self.run_cwl()

        expected_output = {
            'pairs': [
            {
                "pair_id": "Sample10.Sample9", "tumor_id": "Sample10", "normal_id": "Sample9",
                "tmb_maf": OFile(name = "Sample10.Sample9.tmb.maf", size = 1083, hash = "1e5f7f3a6b972c5cd7163d404ecaf99c62311680", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample10.Sample9.tmb.tsv", size = 39, hash = "8d5391d268fb80abaadcfcd158c1d5a206118476", dir = output_dir)
            },
            {
                "pair_id": "Sample12.Sample11", "tumor_id": "Sample12", "normal_id": "Sample11",
                "tmb_maf": OFile(name = "Sample12.Sample11.tmb.maf", size = 3162, hash = "b9764245546db40d1c4138b0a8fd64849b26dbb2", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample12.Sample11.tmb.tsv", size = 39, hash = "0ba713f227e506e35d5af2928fffd15a325d6271", dir = output_dir)
            },
            {
                "pair_id": "Sample14.Sample13", "tumor_id": "Sample14", "normal_id": "Sample13",
                "tmb_maf": OFile(name = "Sample14.Sample13.tmb.maf", size = 1089, hash = "27bced89b50101637da14fec6e204cec0b4bd7f6", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample14.Sample13.tmb.tsv", size = 39, hash = "aa9748f5aa16e01817de2b39b500ac429cc5ca47", dir = output_dir)
            },
            {
                "pair_id": "Sample16.Sample15", "tumor_id": "Sample16", "normal_id": "Sample15",
                "tmb_maf": OFile(name = "Sample16.Sample15.tmb.maf", size = 1473, hash = "4a0487d64c305473e3e00ab475e1f8d403d68caa", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample16.Sample15.tmb.tsv", size = 39, hash = "097d731c93099078d120a1f30cb0e5970a6b003d", dir = output_dir)
            },
            {
                "pair_id": "Sample18.Sample17", "tumor_id": "Sample18", "normal_id": "Sample17",
                "tmb_maf": OFile(name = "Sample18.Sample17.tmb.maf", size = 31897, hash = "d076cb53fae50c6fca195a9999eefc9880192008", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample18.Sample17.tmb.tsv", size = 41, hash = "f3e057e64ac751e8bb314ba5d9d1d5a83a883053", dir = output_dir)
            },
            {
                "pair_id": "Sample20.Sample19", "tumor_id": "Sample20", "normal_id": "Sample19",
                "tmb_maf": OFile(name = "Sample20.Sample19.tmb.maf", size = 1304, hash = "8064bb5f27b6666692b5cecbc4ee36b12fcac9a3", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample20.Sample19.tmb.tsv", size = 39, hash = "87073dc756ead52fff1b5356c8b5cad6f56926d5", dir = output_dir)
            },
            {
                "pair_id": "Sample22.Sample21", "tumor_id": "Sample22", "normal_id": "Sample21",
                "tmb_maf": OFile(name = "Sample22.Sample21.tmb.maf", size = 13946, hash = "8877725be00bc639e8117a462e4f8d1427426778", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample22.Sample21.tmb.tsv", size = 40, hash = "4b3b7bb031164a2a8c5c90bb6ef7378091f1c943", dir = output_dir)
            },
            {
                "pair_id": "Sample24.Sample23", "tumor_id": "Sample24", "normal_id": "Sample23",
                "tmb_maf": OFile(name = "Sample24.Sample23.tmb.maf", size = 3383, hash = "50535e40b1c0c238f17ef2b7759d31a1ed139b77", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample24.Sample23.tmb.tsv", size = 39, hash = "92250bac529328be0d99a14bc43943ddb9da0989", dir = output_dir)
            },
            {
                "pair_id": "Sample26.Sample25", "tumor_id": "Sample26", "normal_id": "Sample25",
                "tmb_maf": OFile(name = "Sample26.Sample25.tmb.maf", size = 6390, hash = "a18318054c6f1121de0e6d7a8895b0c11d9978d4", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample26.Sample25.tmb.tsv", size = 40, hash = "47c5db8fea83bded42ba5b5025a23586189897bf", dir = output_dir)
            },
            {
                "pair_id": "Sample28.Sample27", "tumor_id": "Sample28", "normal_id": "Sample27",
                "tmb_maf": OFile(name = "Sample28.Sample27.tmb.maf", size = 7448, hash = "5b96aaef5d815b4713128123ec6e2b9d51495ea8", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample28.Sample27.tmb.tsv", size = 40, hash = "fbe3ad09f2b9ff72bfcf411a0cf30460219859a6", dir = output_dir)
            },
            {
                "pair_id": "Sample30.Sample29", "tumor_id": "Sample30", "normal_id": "Sample29",
                "tmb_maf": OFile(name = "Sample30.Sample29.tmb.maf", size = 2789, hash = "d1c68680b9118ef234d9692fd5f2def6aa89d763", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample30.Sample29.tmb.tsv", size = 39, hash = "31a3f86b06e206fcdc78aac0d846bf5f7d38d84e", dir = output_dir)
            },
            {
                "pair_id": "Sample32.Sample31", "tumor_id": "Sample32", "normal_id": "Sample31",
                "tmb_maf": OFile(name = "Sample32.Sample31.tmb.maf", size = 1932, hash = "a090868ecdcd90a44563141ff11e60ee3f4c1823", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample32.Sample31.tmb.tsv", size = 39, hash = "887fe7fa431acd1eb5fd8f160fb4d002ed1c1c10", dir = output_dir)
            },
            {
                "pair_id": "Sample34.Sample33", "tumor_id": "Sample34", "normal_id": "Sample33",
                "tmb_maf": OFile(name = "Sample34.Sample33.tmb.maf", size = 2991, hash = "abe7941bc9051d1946c1ae3a60f4b0e72b7e68cf", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample34.Sample33.tmb.tsv", size = 39, hash = "1758446148033024c059a0d288ebe03e0d7a797e", dir = output_dir)
            },
            {
                "pair_id": "Sample50.Sample49", "tumor_id": "Sample50", "normal_id": "Sample49",
                "tmb_maf": OFile(name = "Sample50.Sample49.tmb.maf", size = 18917, hash = "d43ff02fdfbcdad57a140866bf95c0c76dacdeeb", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample50.Sample49.tmb.tsv", size = 40, hash = "429ac67f26f759fd85d930bf18a201d7ad7760af", dir = output_dir)
            },
            {
                "pair_id": "Sample52.Sample51", "tumor_id": "Sample52", "normal_id": "Sample51",
                "tmb_maf": OFile(name = "Sample52.Sample51.tmb.maf", size = 1718, hash = "880020370fa72dbc5bf1ac63ec5fc039b4d28fc7", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample52.Sample51.tmb.tsv", size = 39, hash = "8a2b8e106762ad27669572881a7084aa60fd5323", dir = output_dir)
            },
            {
                "pair_id": "Sample54.Sample53", "tumor_id": "Sample54", "normal_id": "Sample53",
                "tmb_maf": OFile(name = "Sample54.Sample53.tmb.maf", size = 879, hash = "00f0edebec8317ed9f379e518f59590af083ba82", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample54.Sample53.tmb.tsv", size = 36, hash = "76a00ed0315504c795eade8e326d21d06b50b67c", dir = output_dir)
            },
            {
                "pair_id": "Sample62.Sample61", "tumor_id": "Sample62", "normal_id": "Sample61",
                "tmb_maf": OFile(name = "Sample62.Sample61.tmb.maf", size = 1045, hash = "4a41ff19ec6c867df1bbebd0a44babc737c51ba2", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample62.Sample61.tmb.tsv", size = 39, hash = "139810642faa47f70e6d5a357d8fbfc5520cab36", dir = output_dir)
            },
            {
                "pair_id": "Sample70.Sample69", "tumor_id": "Sample70", "normal_id": "Sample69",
                "tmb_maf": OFile(name = "Sample70.Sample69.tmb.maf", size = 17120, hash = "5837551fdd9e95ba8cba85d4150f966e4175aceb", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample70.Sample69.tmb.tsv", size = 40, hash = "c21656cce24af8a10ae3d534e406e705abba2c86", dir = output_dir)
            },
            {
                "pair_id": "Sample74.Sample73", "tumor_id": "Sample74", "normal_id": "Sample73",
                "tmb_maf": OFile(name = "Sample74.Sample73.tmb.maf", size = 8630, hash = "320066d54a202e8abf840dcfe66ebc07774b9b0d", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample74.Sample73.tmb.tsv", size = 40, hash = "c2d30f9498709d2d3c0f5733af1eaaf308d622be", dir = output_dir)
            },
            {
                "pair_id": "Sample80.Sample79", "tumor_id": "Sample80", "normal_id": "Sample79",
                "tmb_maf": OFile(name = "Sample80.Sample79.tmb.maf", size = 3593, hash = "c455da921b4da3e91e8bbe0c8bd06bb0f414debc", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample80.Sample79.tmb.tsv", size = 40, hash = "db19e00007f3fb805373d7dbf5e10539b373fa33", dir = output_dir)
            },
            {
                "pair_id": "Sample84.Sample83", "tumor_id": "Sample84", "normal_id": "Sample83",
                "tmb_maf": OFile(name = "Sample84.Sample83.tmb.maf", size = 4456, hash = "aac9da438d6b909000b36ab44ded93da5818aeec", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample84.Sample83.tmb.tsv", size = 40, hash = "2159332cb06a280cde57f0aa07b011da9e0b3e65", dir = output_dir)
            },
            {
                "pair_id": "Sample90.Sample89", "tumor_id": "Sample90", "normal_id": "Sample89",
                "tmb_maf": OFile(name = "Sample90.Sample89.tmb.maf", size = 1098, hash = "4a8cfb06d53fed560f1e9f05936069f11d57a8e8", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample90.Sample89.tmb.tsv", size = 39, hash = "e4a9595806740ee2224bbcc6c9a34cf398a87e33", dir = output_dir)
            },
            {
                "pair_id": "Sample92.Sample91", "tumor_id": "Sample92", "normal_id": "Sample91",
                "tmb_maf": OFile(name = "Sample92.Sample91.tmb.maf", size = 1088, hash = "7ad7be1538f618ef3462e7f9f8cc7173b462199e", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample92.Sample91.tmb.tsv", size = 39, hash = "0427bbe4c82ddd013d6c2bd489484dd6a6541fce", dir = output_dir)
            },
            ]
            }
        self.assertCWLDictEqual(output_json, expected_output)

        expected_values = {
            'Sample10': '0.8239',
            'Sample12': '9.0627',
            'Sample14': '0.8239',
            'Sample16': '2.4716',
            'Sample18': '120.2864',
            'Sample20': '1.6478',
            'Sample22': '51.0805',
            'Sample24': '9.8866',
            'Sample26': '21.4209',
            'Sample28': '25.5403',
            'Sample30': '7.4149',
            'Sample32': '4.1194',
            'Sample34': '8.2388',
            'Sample50': '70.0297',
            'Sample52': '3.2955',
            'Sample54': '0.0',
            'Sample62': '0.8239',
            'Sample70': '63.4387',
            'Sample74': '30.4835',
            'Sample80': '10.7104',
            'Sample84': '14.0059',
            'Sample90': '0.8239',
            'Sample92': '0.8239'}

        values = self.getAllSampleFileValues(
            filepaths = [ pair["tmb_tsv"]["path"] for pair in expected_output["pairs"] ],
            value_fieldname = "CMO_TMB_SCORE",
            sample_fieldname = "SampleID"
            )

        self.assertDictEqual(values, expected_values)

if __name__ == "__main__":
    unittest.main()
