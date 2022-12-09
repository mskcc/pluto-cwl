#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case for the samples_fillout_index_batch_workflow cwl


$ CWL_ENGINE=Toil PRINT_TESTNAME=T python3 tests/test_samples_fillout_index_batch_workflow_cwl.py TestSamplesFilloutIndexBatch.test_three_groups
"""
import os
import sys
import unittest
from typing import Dict, Tuple

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.cwlFile import CWLFile
from pluto.plutoTestCase import PlutoTestCase
from pluto.plutoPreRunTestCase import PlutoPreRunTestCase
from pluto.settings import DATA_SETS
from pluto.serializer import OFile
sys.path.pop(0)


sample1_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample1.FillOutUnitTest01.muts.maf')
sample2_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample2.FillOutUnitTest01.muts.maf')
sample3_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample3.FillOutUnitTest01.muts.maf')
sample4_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample4.FillOutUnitTest01.muts.maf')
sample5_maf = os.path.join(DATA_SETS['Fillout01']['MAF_DIR'], 'Sample5.FillOutUnitTest01.muts.maf')

sample1_bam = os.path.join(DATA_SETS['Fillout01']['BAM_DIR'], 'Sample1.UnitTest01.bam')
sample2_bam = os.path.join(DATA_SETS['Fillout01']['BAM_DIR'], 'Sample2.UnitTest01.bam')
sample3_bam = os.path.join(DATA_SETS['Fillout01']['BAM_DIR'], 'Sample3.UnitTest01.bam')
sample4_bam = os.path.join(DATA_SETS['Fillout01']['BAM_DIR'], 'Sample4.UnitTest01.bam')
sample5_bam = os.path.join(DATA_SETS['Fillout01']['BAM_DIR'], 'Sample5.UnitTest01.bam')

    # sample_group1 = [
    #         {
    #             "sample_id": "Sample1",
    #             "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
    #             "sample_type": "research",
    #             "prefilter": False,
    #             "maf_file": { "class": "File", "path": self.sample1_maf },
    #             "bam_file": { "class": "File", "path": self.sample1_bam }
    #         },
    #         {
    #             "sample_id": "Sample2",
    #             "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
    #             "sample_type": "research",
    #             "prefilter": False,
    #             "maf_file": { "class": "File", "path": self.sample2_maf },
    #             "bam_file": { "class": "File", "path": self.sample2_bam }
    #         },
    #         {
    #             "sample_id": "Sample3",
    #             "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
    #             "sample_type": "clinical",
    #             "prefilter": False,
    #             "maf_file": { "class": "File", "path": self.sample3_maf },
    #             "bam_file": { "class": "File", "path": self.sample3_bam }
    #         },
    #         {
    #             "sample_id": "Sample4",
    #             "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
    #             "sample_type": "clinical",
    #             "prefilter": False,
    #             "maf_file": { "class": "File", "path": self.sample4_maf },
    #             "bam_file": { "class": "File", "path": self.sample4_bam }
    #         },
    #         {
    #             "sample_id": "Sample5",
    #             "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
    #             "sample_type": "research",
    #             "prefilter": False,
    #             "maf_file": { "class": "File", "path": self.sample5_maf },
    #             "bam_file": { "class": "File", "path": self.sample5_bam }
    #         },
    #     ]


class TestSamplesFilloutIndexBatch1Group(PlutoPreRunTestCase):
    
    cwl_file = CWLFile('samples_fillout_index_batch_workflow.cwl')
    
    def setUp(self):
        super().setUp()
        self.runner_args['use_cache'] = False # do not use cache for samples fillout workflow it breaks on split_vcf_to_mafs

    def setUpRun(self):
        sample_group1 = [
            {
                "sample_id": "Sample1",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "research",
                "prefilter": True,
                "maf_file": { "class": "File", "path": sample1_maf },
                "bam_file": { "class": "File", "path": sample1_bam }
            },
            {
                "sample_id": "Sample2",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "research",
                "prefilter": True,
                "maf_file": { "class": "File", "path": sample2_maf },
                "bam_file": { "class": "File", "path": sample2_bam }
            },
            {
                "sample_id": "Sample3",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "clinical",
                "prefilter": False,
                "maf_file": { "class": "File", "path": sample3_maf },
                "bam_file": { "class": "File", "path": sample3_bam }
            },
        ]

        self.input = {
            "sample_groups": [sample_group1],
            "fillout_output_fname": 'output.maf',
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']},
        }

        output_json, output_dir = self.run_cwl()

        return(output_json, output_dir)
    
    def getExpected(self, output_dir):
        return({
            'output_file': OFile(name = 'output.maf', dir = output_dir),
            'filtered_file': OFile(name = 'output.filtered.maf', dir = output_dir),
            'portal_file': OFile(name = 'data_mutations_extended.txt', dir = output_dir),
            'uncalled_file': OFile(name = 'data_mutations_uncalled.txt', dir = output_dir),
        })
    
    # # # # # # # # # # #
    # # # # # # # # # # #

    def test_CWLDictEqual(self):
        """
        Test case for running the fillout workflow on a number of samples, each with a bam and maf
        """
        # file contents are inconsistent so strip some keys from the output dict
        strip_related_keys = [
        ('basename', 'output.maf', ['size', 'checksum']),
        ('basename', 'output.filtered.maf', ['size', 'checksum']),
        ('basename', 'data_mutations_extended.txt', ['size', 'checksum']),
        ('basename', 'data_mutations_uncalled.txt', ['size', 'checksum'])
        ]
        self.assertCWLDictEqual(
                self.res.output, 
                self.res.expected, 
                related_keys = strip_related_keys)

    def test_output_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['output_file']).path, 147)

    def test_output_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['output_file']).path, "4732e626d2859e4c2e8a7d4eeca0e0f4")
    

    def test_filtered_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['filtered_file']).path, 96)

    def test_filtered_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['filtered_file']).path, "f934e6bd6f1767372b9737d3865e9f0b")

    def test_portal_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['portal_file']).path, 70)

    def test_portal_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['portal_file']).path, "10f4469d0128b6e0bf9e1ef315feb08c")

    def test_uncalled_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['uncalled_file']).path, 26)

    def test_uncalled_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['uncalled_file']).path, "f996e92adc6d1fecb946533a9f23ae99")

    def test_portal_output_path_num_muts(self):
        self.assertEqualNumMutations([
            OFile.init_dict(self.res.output['portal_file']).path, 
            OFile.init_dict(self.res.output['uncalled_file']).path, 
            ], 
            OFile.init_dict(self.res.output['filtered_file']).path)

    def test_output_file_fields(self):
        self.assertMutFieldContains(
            OFile.init_dict(self.res.output['output_file']).path,
            "Tumor_Sample_Barcode", ["Sample1", "Sample2", "Sample3"], containsAll = True)

    def test_portal_output_path_fields(self):
        self.assertMutFieldDoesntContain(
            OFile.init_dict(self.res.output['portal_file']).path,
            "Amino_Acid_Change", [""])

    def test_uncalled_output_path_fields(self):
        self.assertMutFieldDoesntContain(
            OFile.init_dict(self.res.output['uncalled_file']).path,
            "Amino_Acid_Change", [""])









    # TODO: update this test case so it works so that we can check the outputs with/without prefilter and clinical/research sample types
    # NOTE: really need to get framework to run CWL in the 'setup' method so we can check all the assertions in parallel without running it over each time
    # def test_no_prefilter(self):
    #     """
    #     HGVSp_Short Amino_Acid_Change
    #     """
    #     sample_group1 = [
    #         {
    #             "sample_id": "Sample1",
    #             "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
    #             "sample_type": "research",
    #             "prefilter": False,
    #             "maf_file": { "class": "File", "path": self.sample1_maf },
    #             "bam_file": { "class": "File", "path": self.sample1_bam }
    #         },
    #         {
    #             "sample_id": "Sample2",
    #             "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
    #             "sample_type": "research",
    #             "prefilter": False,
    #             "maf_file": { "class": "File", "path": self.sample2_maf },
    #             "bam_file": { "class": "File", "path": self.sample2_bam }
    #         },
    #         {
    #             "sample_id": "Sample3",
    #             "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
    #             "sample_type": "clinical",
    #             "prefilter": False,
    #             "maf_file": { "class": "File", "path": self.sample3_maf },
    #             "bam_file": { "class": "File", "path": self.sample3_bam }
    #         },
    #         {
    #             "sample_id": "Sample4",
    #             "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
    #             "sample_type": "clinical",
    #             "prefilter": False,
    #             "maf_file": { "class": "File", "path": self.sample4_maf },
    #             "bam_file": { "class": "File", "path": self.sample4_bam }
    #         },
    #         {
    #             "sample_id": "Sample5",
    #             "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
    #             "sample_type": "research",
    #             "prefilter": False,
    #             "maf_file": { "class": "File", "path": self.sample5_maf },
    #             "bam_file": { "class": "File", "path": self.sample5_bam }
    #         },
    #     ]

    #     self.input = {
    #         "sample_groups": [sample_group1],
    #         "fillout_output_fname": 'output.maf',
    #         "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']},
    #     }

    #     output_json, output_dir = self.run_cwl()
    #     output_file = os.path.join(output_dir,'output.maf')
    #     filtered_output_path = os.path.join(output_dir,'output.filtered.maf')
    #     portal_output_path = os.path.join(output_dir,'data_mutations_extended.txt')
    #     uncalled_output_path = os.path.join(output_dir,'data_mutations_uncalled.txt')

    #     expected_output = {
    #         'output_file': OFile(name = 'output.maf', dir = output_dir),
    #         'filtered_file': OFile(name = 'output.filtered.maf', dir = output_dir),
    #         'portal_file': OFile(name = 'data_mutations_extended.txt', dir = output_dir),
    #         'uncalled_file': OFile(name = 'data_mutations_uncalled.txt', dir = output_dir),
    #     }

    #     # file contents are inconsistent so strip some keys from the output dict
    #     strip_related_keys = [
    #     ('basename', 'output.maf', ['size', 'checksum']),
    #     ('basename', 'output.filtered.maf', ['size', 'checksum']),
    #     ('basename', 'data_mutations_extended.txt', ['size', 'checksum']),
    #     ('basename', 'data_mutations_uncalled.txt', ['size', 'checksum'])
    #     ]
    #     self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)

    #     self.assertNumMutationsHash(output_file, 310, '18fafe6dd335cb62f515e0323e6b74b2')
    #     self.assertNumMutationsHash(filtered_output_path, 225, '450b97a2b93ed9421c141837f99240ce')
    #     self.assertNumMutationsHash(portal_output_path, 159, '82c2ab2962f782494ddd87886f1ff03b')
    #     self.assertNumMutationsHash(uncalled_output_path, 66, 'e866fffc38c7d0b2602617973e496039')
    #     self.assertEqualNumMutations([portal_output_path, uncalled_output_path], filtered_output_path)
    #     self.assertMutFieldContains(output_file, "Tumor_Sample_Barcode", ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5"], containsAll = True)
    #     self.assertMutFieldDoesntContain(portal_output_path, "Amino_Acid_Change", [""])
    #     self.assertMutFieldDoesntContain(uncalled_output_path, "Amino_Acid_Change", [""])










class Dummy(object):
    def test_two_groups(self):
        sample_group1 = [
            {
                "sample_id": "Sample1",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "research",
                "prefilter": True,
                "maf_file": { "class": "File", "path": self.sample1_maf },
                "bam_file": { "class": "File", "path": self.sample1_bam }
            },
            {
                "sample_id": "Sample2",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "research",
                "prefilter": True,
                "maf_file": { "class": "File", "path": self.sample2_maf },
                "bam_file": { "class": "File", "path": self.sample2_bam }
            },
            {
                "sample_id": "Sample3",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "clinical",
                "prefilter": False,
                "maf_file": { "class": "File", "path": self.sample3_maf },
                "bam_file": { "class": "File", "path": self.sample3_bam }
            },
        ]

        sample_group2 = [
            {
                "sample_id": "Sample4",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "clinical",
                "prefilter": False,
                "maf_file": { "class": "File", "path": self.sample4_maf },
                "bam_file": { "class": "File", "path": self.sample4_bam }
            },
            {
                "sample_id": "Sample5",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "research",
                "prefilter": True,
                "maf_file": { "class": "File", "path": self.sample5_maf },
                "bam_file": { "class": "File", "path": self.sample5_bam }
            },
        ]

        self.input = {
            "sample_groups": [sample_group1, sample_group2],
            "fillout_output_fname": 'output.maf',
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']},
        }

        output_json, output_dir = self.run_cwl()
        output_file = os.path.join(output_dir,'output.maf')
        filtered_output_path = os.path.join(output_dir,'output.filtered.maf')
        portal_output_path = os.path.join(output_dir,'data_mutations_extended.txt')
        uncalled_output_path = os.path.join(output_dir,'data_mutations_uncalled.txt')

        expected_output = {
            'output_file': OFile(name = 'output.maf', dir = output_dir),
            'filtered_file': OFile(name = 'output.filtered.maf', dir = output_dir),
            'portal_file': OFile(name = 'data_mutations_extended.txt', dir = output_dir),
            'uncalled_file': OFile(name = 'data_mutations_uncalled.txt', dir = output_dir),
        }

        # file contents are inconsistent so strip some keys from the output dict
        strip_related_keys = [
        ('basename', 'output.maf', ['size', 'checksum']),
        ('basename', 'output.filtered.maf', ['size', 'checksum']),
        ('basename', 'data_mutations_extended.txt', ['size', 'checksum']),
        ('basename', 'data_mutations_uncalled.txt', ['size', 'checksum'])
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)

        self.assertNumMutationsHash(output_file, 235, '4e4c91ef129a853a35b86f7fa6f1268a')
        self.assertNumMutationsHash(filtered_output_path, 150, '8397a12302977db14e798a1b2e3ba151')
        self.assertNumMutationsHash(portal_output_path, 120, 'dc0941e362f97df67508cd70d2d8a76e')
        self.assertNumMutationsHash(uncalled_output_path, 30, '3b78a67aaa0429cc577578910faf5c10')
        self.assertEqualNumMutations([portal_output_path, uncalled_output_path], filtered_output_path)
        self.assertMutFieldContains(output_file, "Tumor_Sample_Barcode", ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5"], containsAll = True)




    def test_three_groups(self):
        """
        NOTE: This one includes Singleton case! most important
        Three groups, one of which contains a single sample (singleton)
        """
        sample_group1 = [
            {
                "sample_id": "Sample1",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "research",
                "prefilter": True,
                "maf_file": { "class": "File", "path": self.sample1_maf },
                "bam_file": { "class": "File", "path": self.sample1_bam }
            },
            {
                "sample_id": "Sample2",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "research",
                "prefilter": True,
                "maf_file": { "class": "File", "path": self.sample2_maf },
                "bam_file": { "class": "File", "path": self.sample2_bam }
            },
        ]

        sample_group2 = [
            {
                "sample_id": "Sample3",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "clinical",
                "prefilter": False,
                "maf_file": { "class": "File", "path": self.sample3_maf },
                "bam_file": { "class": "File", "path": self.sample3_bam }
            },
            {
                "sample_id": "Sample4",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "clinical",
                "prefilter": False,
                "maf_file": { "class": "File", "path": self.sample4_maf },
                "bam_file": { "class": "File", "path": self.sample4_bam }
            }
        ]

        # Singleton sample; no DMP clinical matches
        sample_group3 = [
            {
                "sample_id": "Sample5",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "research",
                "prefilter": True,
                "maf_file": { "class": "File", "path": self.sample5_maf },
                "bam_file": { "class": "File", "path": self.sample5_bam }
            },
        ]

        self.input = {
            "sample_groups": [sample_group1, sample_group2, sample_group3],
            "fillout_output_fname": 'output.maf',
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']},
        }

        output_json, output_dir = self.run_cwl()
        output_file = os.path.join(output_dir,'output.maf')
        filtered_output_path = os.path.join(output_dir,'output.filtered.maf')
        portal_output_path = os.path.join(output_dir,'data_mutations_extended.txt')
        uncalled_output_path = os.path.join(output_dir,'data_mutations_uncalled.txt')

        expected_output = {
            'output_file': OFile(name = 'output.maf', dir = output_dir),
            'filtered_file': OFile(name = 'output.filtered.maf', dir = output_dir),
            'portal_file': OFile(name = 'data_mutations_extended.txt', dir = output_dir),
            'uncalled_file': OFile(name = 'data_mutations_uncalled.txt', dir = output_dir),
        }

        # file contents are inconsistent so strip some keys from the output dict
        strip_related_keys = [
        ('basename', 'output.maf', ['size', 'checksum']),
        ('basename', 'output.filtered.maf', ['size', 'checksum']),
        ('basename', 'data_mutations_extended.txt', ['size', 'checksum']),
        ('basename', 'data_mutations_uncalled.txt', ['size', 'checksum'])
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)

        self.assertNumMutationsHash(output_file, 188, '2f4c5e5cb13430f456bbc41a0a93dc41')
        self.assertNumMutationsHash(filtered_output_path, 126, '3dda4952d2ae396079155b4bc8cc276f')
        self.assertNumMutationsHash(portal_output_path, 108, 'f4b365c0e6be4b975a1443905860ffb6')
        self.assertNumMutationsHash(uncalled_output_path, 18, '7405745e675ea524d8e16e09b9bb749d')
        self.assertEqualNumMutations([portal_output_path, uncalled_output_path], filtered_output_path)
        self.assertMutFieldContains(output_file, "Tumor_Sample_Barcode", ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5"], containsAll = True)






    def test_four_groups(self):
        """
        Four groups, two of which contains a single sample (singleton)
        """
        sample_group1 = [
            {
                "sample_id": "Sample1",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "research",
                "prefilter": True,
                "maf_file": { "class": "File", "path": self.sample1_maf },
                "bam_file": { "class": "File", "path": self.sample1_bam }
            },
            {
                "sample_id": "Sample2",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "research",
                "prefilter": True,
                "maf_file": { "class": "File", "path": self.sample2_maf },
                "bam_file": { "class": "File", "path": self.sample2_bam }
            },
        ]

        sample_group2 = [
            {
                "sample_id": "Sample3",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "clinical",
                "prefilter": False,
                "maf_file": { "class": "File", "path": self.sample3_maf },
                "bam_file": { "class": "File", "path": self.sample3_bam }
            },
        ]

        sample_group3 = [
            {
                "sample_id": "Sample5",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "research",
                "prefilter": True,
                "maf_file": { "class": "File", "path": self.sample5_maf },
                "bam_file": { "class": "File", "path": self.sample5_bam }
            },
        ]

        sample_group4 = [
            {
                "sample_id": "Sample4",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "clinical",
                "prefilter": False,
                "maf_file": { "class": "File", "path": self.sample4_maf },
                "bam_file": { "class": "File", "path": self.sample4_bam }
            }
        ]

        self.input = {
            "sample_groups": [sample_group1, sample_group2, sample_group3, sample_group4],
            "fillout_output_fname": 'output.maf',
            "ref_fasta": {"class": "File", "path": self.DATA_SETS['Proj_08390_G']['REF_FASTA']},
        }

        output_json, output_dir = self.run_cwl()
        output_file = os.path.join(output_dir,'output.maf')
        filtered_output_path = os.path.join(output_dir,'output.filtered.maf')
        portal_output_path = os.path.join(output_dir,'data_mutations_extended.txt')
        uncalled_output_path = os.path.join(output_dir,'data_mutations_uncalled.txt')

        expected_output = {
            'output_file': OFile(name = 'output.maf', dir = output_dir),
            'filtered_file': OFile(name = 'output.filtered.maf', dir = output_dir),
            'portal_file': OFile(name = 'data_mutations_extended.txt', dir = output_dir),
            'uncalled_file': OFile(name = 'data_mutations_uncalled.txt', dir = output_dir),
        }

        # file contents are inconsistent so strip some keys from the output dict
        strip_related_keys = [
        ('basename', 'output.maf', ['size', 'checksum']),
        ('basename', 'output.filtered.maf', ['size', 'checksum']),
        ('basename', 'data_mutations_extended.txt', ['size', 'checksum']),
        ('basename', 'data_mutations_uncalled.txt', ['size', 'checksum'])
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)

        self.assertNumMutationsHash(output_file, 157, 'ed7dcce977a13f360463e45f5a07154b') # , _print = True
        self.assertNumMutationsHash(filtered_output_path, 36, '5ea9c4b66287a100fc90e05619d52364')
        self.assertNumMutationsHash(portal_output_path, 36, '2a9b7ca9c0942c4daf3273c63858db3c')
        self.assertNumMutationsHash(uncalled_output_path, 0, 'd751713988987e9331980363e24189ce')
        self.assertEqualNumMutations([portal_output_path, uncalled_output_path], filtered_output_path)
        self.assertMutFieldContains(output_file, "Tumor_Sample_Barcode", ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5"], containsAll = True)

if __name__ == "__main__":
    unittest.main()
