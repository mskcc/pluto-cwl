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
from pluto.tools import PlutoTestCase, CWLFile, TableReader
from pluto.serializer import OFile
sys.path.pop(0)

# handle for errors arising from python3 -m unittest ...
try:
    import fixtures_tmb as fxt
except ModuleNotFoundError:
    sys.path.insert(0, THIS_DIR)
    import fixtures_tmb as fxt
    sys.path.pop(0)


class TestTmbWorkflow(PlutoTestCase):
    cwl_file = CWLFile('tmb_workflow.cwl')

    def setUp(self):
        super().setUp() # initialize the tmpdir

        # PASS: 6
        # FAIL: 4
        self.maf_lines1 = self.dicts2lines(dict_list = fxt.rows1, comment_list = fxt.comments)

        # PASS: 5
        # FAIL: 5
        self.maf_lines2 = self.dicts2lines(dict_list = fxt.rows2, comment_list = fxt.comments)
        self.maf1 = self.write_table(self.tmpdir, filename = "input1.maf", lines = self.maf_lines1)
        self.maf2 = self.write_table(self.tmpdir, filename = "input2.maf", lines = self.maf_lines2)

    def test_tmb_workflow1(self):
        """
        Test case for running the TMB workflow on multiple files
        """
        self.input = {
            "assay_coverage":  '1000',
            "pairs": [
                {
                    "pair_maf": {
                        "path": self.maf1,
                        "class": "File"
                    },
                    "pair_id": "Sample1-T.Sample1-N",
                    "tumor_id": "Sample1-T",
                    "normal_id": "Sample1-N"
                },
                {
                    "pair_maf": {
                        "path": self.maf2,
                        "class": "File"
                    },
                    "pair_id": "Sample2-T.Sample2-N",
                    "tumor_id": "Sample2-T",
                    "normal_id": "Sample2-N"
                }
                ]
            }

        output_json, output_dir = self.run_cwl()

        expected_output = {
        "pairs":[
            {
                "pair_id": "Sample1-T.Sample1-N",
                "tumor_id": "Sample1-T",
                "normal_id": "Sample1-N",
                "tmb_maf": OFile(name = "Sample1-T.Sample1-N.tmb.maf", size = 352, hash = "b019b5b3c6aba861371c135fe47520c969fab5ae", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample1-T.Sample1-N.tmb.tsv", size = 40, hash = "d6a57cfb5e3001697875e5b5bfae206e0f7f2310", dir = output_dir)
            },
            {
                "pair_id": "Sample2-T.Sample2-N",
                "tumor_id": "Sample2-T",
                "normal_id": "Sample2-N",
                "tmb_maf": OFile(name = "Sample2-T.Sample2-N.tmb.maf", size = 308, hash = "e6cb790e887d606e2762dd73e73fec28ecfff22b", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample2-T.Sample2-N.tmb.tsv", size = 40, hash = "851b303d4237654de54fdfdc60d41ab996f4380d", dir = output_dir)
            }
            ]
            }
        self.assertCWLDictEqual(output_json, expected_output)

        lines = self.read_table(expected_output["pairs"][0]["tmb_tsv"]["path"])
        expected_lines = [
            ['CMO_TMB_SCORE', 'SampleID'],
            ['7000.0', 'Sample1-T']
        ]
        self.assertEqual(lines, expected_lines)

        lines = self.read_table(expected_output["pairs"][1]["tmb_tsv"]["path"])
        expected_lines = [
            ['CMO_TMB_SCORE', 'SampleID'],
            ['6000.0', 'Sample2-T']
        ]
        self.assertEqual(lines, expected_lines)


    def test_tmb_workflow1_2(self):
        # one pair uses a Pooled Normal
        self.input = {
            "assay_coverage":  '1000',
            "pairs": [
                {
                    "pair_maf": {
                        "path": self.maf1,
                        "class": "File"
                    },
                    "pair_id": "Sample1-T.Sample1-PooledNormal",
                    "tumor_id": "Sample1-T",
                    "normal_id": "Sample1-PooledNormal"
                },
                {
                    "pair_maf": {
                        "path": self.maf2,
                        "class": "File"
                    },
                    "pair_id": "Sample2-T.Sample2-N",
                    "tumor_id": "Sample2-T",
                    "normal_id": "Sample2-N"
                }
                ]
            }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            "pairs":[
                {
                    "pair_id": "Sample1-T.Sample1-PooledNormal",
                    "tumor_id": "Sample1-T",
                    "normal_id": "Sample1-PooledNormal",
                    "tmb_maf": OFile(name = "Sample1-T.Sample1-PooledNormal.tmb.maf", size = 352, hash = "b019b5b3c6aba861371c135fe47520c969fab5ae", dir = output_dir),
                    "tmb_tsv": OFile(name = "Sample1-T.Sample1-PooledNormal.tmb.tsv", size = 36, hash = "84ba7be1f81eba42d9845ccb18cae9d4b3d6a21f", dir = output_dir)
                },
                {
                    "pair_id": "Sample2-T.Sample2-N",
                    "tumor_id": "Sample2-T",
                    "normal_id": "Sample2-N",
                    "tmb_maf": OFile(name = "Sample2-T.Sample2-N.tmb.maf", size = 308, hash = "e6cb790e887d606e2762dd73e73fec28ecfff22b", dir = output_dir),
                    "tmb_tsv": OFile(name = "Sample2-T.Sample2-N.tmb.tsv", size = 40, hash = "851b303d4237654de54fdfdc60d41ab996f4380d", dir = output_dir)
                }
                ]
            }
        self.assertCWLDictEqual(output_json, expected_output)

        # !!!!!!!!
        # Sample1 gets value NA because its paired with a pooled normal !!
        lines = self.read_table(expected_output["pairs"][0]["tmb_tsv"]["path"])
        expected_lines = [
            ['CMO_TMB_SCORE', 'SampleID'],
            ['NA', 'Sample1-T']
        ]
        self.assertEqual(lines, expected_lines)
        # !!!!!!!!

        lines = self.read_table(expected_output["pairs"][1]["tmb_tsv"]["path"])
        expected_lines = [
            ['CMO_TMB_SCORE', 'SampleID'],
            ['6000.0', 'Sample2-T']
        ]
        self.assertEqual(lines, expected_lines)

    def test_tmb_workflow2(self):
        """
        Test case for using a single input pair maf
        """
        self.input = {
            "assay_coverage":  '1000',
            "pairs": [
                {
                    "pair_maf": {"path": self.maf1, "class": "File"},
                    "pair_id": "Sample1-T.Sample1-N",
                    "tumor_id": "Sample1-T",
                    "normal_id": "Sample1-N"
                }
                ]
            }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'pairs': [{
                "pair_id": "Sample1-T.Sample1-N",
                "tumor_id": "Sample1-T",
                "normal_id": "Sample1-N",
                "tmb_maf": OFile(name = "Sample1-T.Sample1-N.tmb.maf", size = 352, hash = "b019b5b3c6aba861371c135fe47520c969fab5ae", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample1-T.Sample1-N.tmb.tsv", size = 40, hash = "d6a57cfb5e3001697875e5b5bfae206e0f7f2310", dir = output_dir)
            }]
            }
        self.assertCWLDictEqual(output_json, expected_output)

        lines = self.read_table(expected_output["pairs"][0]["tmb_tsv"]["path"])
        expected_lines = [
            ['CMO_TMB_SCORE', 'SampleID'],
            ['7000.0', 'Sample1-T']
        ]
        self.assertEqual(lines, expected_lines)


    def test_tmb_workflow3(self):
        """
        Test case with a single real maf file
        """
        self.input = {
            "assay_coverage":  '1000',
            "pairs": [
                {
                    "pair_maf": {
                        "path": os.path.join(self.DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf"),
                        "class": "File"
                    },
                    "pair_id": "Sample1-T.Sample1-N",
                    "tumor_id": "Sample1-T",
                    "normal_id": "Sample1-N"
                }
                ]
            }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'pairs': [{
                "pair_id": "Sample1-T.Sample1-N",
                "tumor_id": "Sample1-T",
                "normal_id": "Sample1-N",
                "tmb_maf": OFile(name = "Sample1-T.Sample1-N.tmb.maf", size = 519440, hash = "809c3c1ac3bb750aebf22ee2f95a5ebafd41e98f", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample1-T.Sample1-N.tmb.tsv", size = 42, hash = "8156b9e7a0602ddd7710f002ef9385237a82c5d0", dir = output_dir)
            }]
            }
        self.assertCWLDictEqual(output_json, expected_output)

        lines = self.read_table(expected_output["pairs"][0]["tmb_tsv"]["path"])
        expected_lines = [
            ['CMO_TMB_SCORE', 'SampleID'],
            ['475000.0', 'Sample1-T']
        ]
        self.assertEqual(lines, expected_lines)

    def test_tmb_workflow4(self):
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
                "tmb_maf": OFile(name = "Sample10.Sample9.tmb.maf", size = 435695, hash = "bf956cd3c611a398e4773061074e34ba7a19de0e", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample10.Sample9.tmb.tsv", size = 41, hash = "88c366b5c6fb06f042d7ed93dbb11c1b660f8a09", dir = output_dir)
            },
            {
                "pair_id": "Sample12.Sample11", "tumor_id": "Sample12", "normal_id": "Sample11",
                "tmb_maf": OFile(name = "Sample12.Sample11.tmb.maf", size = 349189, hash = "e0ffdde8c70d3b2f9d9b25dee279524d056d7f11", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample12.Sample11.tmb.tsv", size = 41, hash = "b2722238425a666c5251b526719b9f18a802b870", dir = output_dir)
            },
            {
                "pair_id": "Sample14.Sample13", "tumor_id": "Sample14", "normal_id": "Sample13",
                "tmb_maf": OFile(name = "Sample14.Sample13.tmb.maf", size = 175012, hash = "ba440532988647cb861876b55e8f7c7818932924", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample14.Sample13.tmb.tsv", size = 41, hash = "24268f13b478d05fdc42aab378619ac8166a85be", dir = output_dir)
            },
            {
                "pair_id": "Sample16.Sample15", "tumor_id": "Sample16", "normal_id": "Sample15",
                "tmb_maf": OFile(name = "Sample16.Sample15.tmb.maf", size = 132771, hash = "8567df193d795fd228bd4c03f168fda7c3cc61a6", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample16.Sample15.tmb.tsv", size = 41, hash = "7a10e353bd4b85197a51b1a01d97b4d5fde4a0d6", dir = output_dir)
            },
            {
                "pair_id": "Sample18.Sample17", "tumor_id": "Sample18", "normal_id": "Sample17",
                "tmb_maf": OFile(name = "Sample18.Sample17.tmb.maf", size = 1198414, hash = "f23cba80210c7053de7d7c5419366125e83643a5", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample18.Sample17.tmb.tsv", size = 41, hash = "df324eb65e39e98ba7dfbe708311f6483f2b7a89", dir = output_dir)
            },
            {
                "pair_id": "Sample20.Sample19", "tumor_id": "Sample20", "normal_id": "Sample19",
                "tmb_maf": OFile(name = "Sample20.Sample19.tmb.maf", size = 452313, hash = "6715dfc2887a5315dd74405fe416c77cbdda7589", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample20.Sample19.tmb.tsv", size = 41, hash = "65285f51785040f8b25be89602611466619b8ca6", dir = output_dir)
            },
            {
                "pair_id": "Sample22.Sample21", "tumor_id": "Sample22", "normal_id": "Sample21",
                "tmb_maf": OFile(name = "Sample22.Sample21.tmb.maf", size = 657290, hash = "c5a47f421bdb260e00ff5ec4224a4854cf083f8f", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample22.Sample21.tmb.tsv", size = 41, hash = "ea02f25e8803d6797a07032c6ba506326bc5259a", dir = output_dir)
            },
            {
                "pair_id": "Sample24.Sample23", "tumor_id": "Sample24", "normal_id": "Sample23",
                "tmb_maf": OFile(name = "Sample24.Sample23.tmb.maf", size = 106895, hash = "92c59891e947f241ad1dad88b47b3da2e243c3a6", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample24.Sample23.tmb.tsv", size = 40, hash = "7bf9f34e1b52d03e702708cda81c5236033421a0", dir = output_dir)
            },
            {
                "pair_id": "Sample26.Sample25", "tumor_id": "Sample26", "normal_id": "Sample25",
                "tmb_maf": OFile(name = "Sample26.Sample25.tmb.maf", size = 378363, hash = "af84468dfc91a50fc2fafb84b060128d022fc2a6", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample26.Sample25.tmb.tsv", size = 41, hash = "05ea960219516542703d2a7254306be241ab8fc3", dir = output_dir)
            },
            {
                "pair_id": "Sample28.Sample27", "tumor_id": "Sample28", "normal_id": "Sample27",
                "tmb_maf": OFile(name = "Sample28.Sample27.tmb.maf", size = 704143, hash = "7ef7cdff11c12191d1948ac19cfd015b6a258ebb", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample28.Sample27.tmb.tsv", size = 41, hash = "41b7b9b7af35791e4ba1ba922c28319ff166dd9c", dir = output_dir)
            },
            {
                "pair_id": "Sample30.Sample29", "tumor_id": "Sample30", "normal_id": "Sample29",
                "tmb_maf": OFile(name = "Sample30.Sample29.tmb.maf", size = 144459, hash = "f9a783df443cfc67dca7f9b110c373b787029128", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample30.Sample29.tmb.tsv", size = 40, hash = "02094abb69f268710fdd0b3a413db89824c8c584", dir = output_dir)
            },
            {
                "pair_id": "Sample32.Sample31", "tumor_id": "Sample32", "normal_id": "Sample31",
                "tmb_maf": OFile(name = "Sample32.Sample31.tmb.maf", size = 325388, hash = "bf9107db224d12f3e033da6c52ed8f42384abcc3", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample32.Sample31.tmb.tsv", size = 41, hash = "6f53a8095503affa1ea4acee2ddade5184cd96c9", dir = output_dir)
            },
            {
                "pair_id": "Sample34.Sample33", "tumor_id": "Sample34", "normal_id": "Sample33",
                "tmb_maf": OFile(name = "Sample34.Sample33.tmb.maf", size = 169118, hash = "d4d8704fe1273a2f3973168334531e62760a038a", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample34.Sample33.tmb.tsv", size = 41, hash = "a1251999d4c785af60cd2696a450f64508c8d47c", dir = output_dir)
            },
            {
                "pair_id": "Sample50.Sample49", "tumor_id": "Sample50", "normal_id": "Sample49",
                "tmb_maf": OFile(name = "Sample50.Sample49.tmb.maf", size = 948931, hash = "e9cc623e4838891bf6b6dd477a1544bc2963f6e4", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample50.Sample49.tmb.tsv", size = 41, hash = "24c941afe6749191910a7f03c590edd4c5d1804a", dir = output_dir)
            },
            {
                "pair_id": "Sample52.Sample51", "tumor_id": "Sample52", "normal_id": "Sample51",
                "tmb_maf": OFile(name = "Sample52.Sample51.tmb.maf", size = 1373542, hash = "e42bd19f856662e9b8f613ce32e8b08c59a21da7", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample52.Sample51.tmb.tsv", size = 41, hash = "02ec3d615c357e0e42697d6bbd2ef8f845af7778", dir = output_dir)
            },
            {
                "pair_id": "Sample54.Sample53", "tumor_id": "Sample54", "normal_id": "Sample53",
                "tmb_maf": OFile(name = "Sample54.Sample53.tmb.maf", size = 203312, hash = "4e013c35588c4373bcea2772f5720d341c8f5617", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample54.Sample53.tmb.tsv", size = 41, hash = "9c86431a33c224b75ea31345056e4fd09b9e6c39", dir = output_dir)
            },
            {
                "pair_id": "Sample62.Sample61", "tumor_id": "Sample62", "normal_id": "Sample61",
                "tmb_maf": OFile(name = "Sample62.Sample61.tmb.maf", size = 280385, hash = "2bcb5509921551fc3e403e666343450038a2acf5", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample62.Sample61.tmb.tsv", size = 41, hash = "cde91a97789f4248edd428eb7aa4f0519b1483d0", dir = output_dir)
            },
            {
                "pair_id": "Sample70.Sample69", "tumor_id": "Sample70", "normal_id": "Sample69",
                "tmb_maf": OFile(name = "Sample70.Sample69.tmb.maf", size = 1591529, hash = "24e7b6c07a6da7e6bc31e1f77f7514050a6c34c0", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample70.Sample69.tmb.tsv", size = 42, hash = "23c1715e2d10d358262b2fb8a7949bfd3a6b8984", dir = output_dir)
            },
            {
                "pair_id": "Sample74.Sample73", "tumor_id": "Sample74", "normal_id": "Sample73",
                "tmb_maf": OFile(name = "Sample74.Sample73.tmb.maf", size = 634232, hash = "8ea2c386032a01a1e0a11631cd002ac88037caa0", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample74.Sample73.tmb.tsv", size = 40, hash = "b7b7fb3b96b5bd02dadc9998f896320379ea19db", dir = output_dir)
            },
            {
                "pair_id": "Sample80.Sample79", "tumor_id": "Sample80", "normal_id": "Sample79",
                "tmb_maf": OFile(name = "Sample80.Sample79.tmb.maf", size = 242884, hash = "6b12f3ec77d76d3d15ec2d837e6842cff9ca83f6", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample80.Sample79.tmb.tsv", size = 41, hash = "5ccc6c6a1c658e1b712fa612540239aed77058c6", dir = output_dir)
            },
            {
                "pair_id": "Sample84.Sample83", "tumor_id": "Sample84", "normal_id": "Sample83",
                "tmb_maf": OFile(name = "Sample84.Sample83.tmb.maf", size = 428091, hash = "da7c20bb6ec5a57640d5c506243d2a867769a3d1", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample84.Sample83.tmb.tsv", size = 41, hash = "15e4d7786a0387ddee263459a6621b9075a8c35a", dir = output_dir)
            },
            {
                "pair_id": "Sample90.Sample89", "tumor_id": "Sample90", "normal_id": "Sample89",
                "tmb_maf": OFile(name = "Sample90.Sample89.tmb.maf", size = 923130, hash = "172f4903c381e511f4d2008b309837308dbcbb14", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample90.Sample89.tmb.tsv", size = 41, hash = "47fb8577a6b01aee6f486c1c1cc1ad4761c2b700", dir = output_dir)
            },
            {
                "pair_id": "Sample92.Sample91", "tumor_id": "Sample92", "normal_id": "Sample91",
                "tmb_maf": OFile(name = "Sample92.Sample91.tmb.maf", size = 242115, hash = "26262d0012d3bed3ba13f9618931aca70030752c", dir = output_dir),
                "tmb_tsv": OFile(name = "Sample92.Sample91.tmb.tsv", size = 41, hash = "c6e022e3fa5335f352be6840da5c40cc91aef0b7", dir = output_dir)
            },
            ]
            }
        self.assertCWLDictEqual(output_json, expected_output)

        expected_values = {
            'Sample10': '312.2503',
            'Sample12': '247.9877',
            'Sample14': '125.2297',
            'Sample16': '101.3372',
            'Sample18': '893.0852',
            'Sample20': '325.4323',
            'Sample22': '506.6858',
            'Sample24': '69.2059',
            'Sample26': '270.2324',
            'Sample28': '519.8678',
            'Sample30': '94.7461',
            'Sample32': '248.8116',
            'Sample34': '106.2804',
            'Sample50': '702.7691',
            'Sample52': '954.8761',
            'Sample54': '141.7072',
            'Sample62': '217.5041',
            'Sample70': '1108.1177',
            'Sample74': '457.253',
            'Sample80': '168.0714',
            'Sample84': '302.3637',
            'Sample90': '688.7631',
            'Sample92': '163.1281'}

        values = self.getAllSampleFileValues(
            filepaths = [ pair["tmb_tsv"]["path"] for pair in expected_output["pairs"] ],
            value_fieldname = "CMO_TMB_SCORE",
            sample_fieldname = "SampleID"
            )

        self.assertDictEqual(values, expected_values)

if __name__ == "__main__":
    unittest.main()
