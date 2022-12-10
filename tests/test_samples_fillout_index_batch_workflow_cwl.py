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


class TestSamplesFilloutIndexBatch2Group(PlutoPreRunTestCase):
    
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
            }]
        sample_group2 = [
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
            "sample_groups": [sample_group1, sample_group2],
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
        self.assertNumMutations(OFile.init_dict(self.res.output['output_file']).path, 68)

    def test_output_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['output_file']).path, "89b2574c5a8ae02ac44e0bbe897a4bf3")
    

    def test_filtered_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['filtered_file']).path, 68)

    def test_filtered_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['filtered_file']).path, "89b2574c5a8ae02ac44e0bbe897a4bf3")

    def test_portal_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['portal_file']).path, 68)

    def test_portal_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['portal_file']).path, "615c2bb7c266ae2c2c66f91daf745983")

    def test_uncalled_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['uncalled_file']).path, 0)

    def test_uncalled_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['uncalled_file']).path, "d751713988987e9331980363e24189ce")

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


class TestSamplesFilloutIndexBatch2Group2(PlutoPreRunTestCase):
    
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

        sample_group2 = [
            {
                "sample_id": "Sample4",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "clinical",
                "prefilter": False,
                "maf_file": { "class": "File", "path": sample4_maf },
                "bam_file": { "class": "File", "path": sample4_bam }
            },
            {
                "sample_id": "Sample5",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "research",
                "prefilter": True,
                "maf_file": { "class": "File", "path": sample5_maf },
                "bam_file": { "class": "File", "path": sample5_bam }
            },
        ]

        self.input = {
            "sample_groups": [sample_group1, sample_group2],
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
        self.assertNumMutations(OFile.init_dict(self.res.output['output_file']).path, 235)

    def test_output_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['output_file']).path, "4e4c91ef129a853a35b86f7fa6f1268a")    

    def test_filtered_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['filtered_file']).path, 150)

    def test_filtered_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['filtered_file']).path, "8397a12302977db14e798a1b2e3ba151")

    def test_portal_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['portal_file']).path, 120)

    def test_portal_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['portal_file']).path, "dc0941e362f97df67508cd70d2d8a76e")

    def test_uncalled_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['uncalled_file']).path, 30)

    def test_uncalled_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['uncalled_file']).path, "3b78a67aaa0429cc577578910faf5c10")

    def test_portal_output_path_num_muts(self):
        self.assertEqualNumMutations([
            OFile.init_dict(self.res.output['portal_file']).path, 
            OFile.init_dict(self.res.output['uncalled_file']).path, 
            ], 
            OFile.init_dict(self.res.output['filtered_file']).path)

    def test_output_file_fields(self):
        self.assertMutFieldContains(
            OFile.init_dict(self.res.output['output_file']).path,
            "Tumor_Sample_Barcode", ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5"], containsAll = True)

    def test_portal_output_path_fields(self):
        self.assertMutFieldDoesntContain(
            OFile.init_dict(self.res.output['portal_file']).path,
            "Amino_Acid_Change", [""])

    def test_uncalled_output_path_fields(self):
        self.assertMutFieldDoesntContain(
            OFile.init_dict(self.res.output['uncalled_file']).path,
            "Amino_Acid_Change", [""])


class TestSamplesFilloutIndexBatch3Group(PlutoPreRunTestCase):
    
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
        ]

        sample_group2 = [
            {
                "sample_id": "Sample3",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "clinical",
                "prefilter": False,
                "maf_file": { "class": "File", "path": sample3_maf },
                "bam_file": { "class": "File", "path": sample3_bam }
            },
            {
                "sample_id": "Sample4",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "clinical",
                "prefilter": False,
                "maf_file": { "class": "File", "path": sample4_maf },
                "bam_file": { "class": "File", "path": sample4_bam }
            }
        ]

        # Singleton sample; no DMP clinical matches
        sample_group3 = [
            {
                "sample_id": "Sample5",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "research",
                "prefilter": True,
                "maf_file": { "class": "File", "path": sample5_maf },
                "bam_file": { "class": "File", "path": sample5_bam }
            },
        ]

        self.input = {
            "sample_groups": [sample_group1, sample_group2, sample_group3],
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
        self.assertNumMutations(OFile.init_dict(self.res.output['output_file']).path, 144)

    def test_output_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['output_file']).path, "85abef967b1d43112da6a026e80f5cea")    

    def test_filtered_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['filtered_file']).path, 144)

    def test_filtered_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['filtered_file']).path, "6614e73fc609a8cfaba4d0a99f9773ea")

    def test_portal_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['portal_file']).path, 126)

    def test_portal_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['portal_file']).path, "683a9892a9727960da82532328b43ae9")

    def test_uncalled_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['uncalled_file']).path, 18)

    def test_uncalled_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['uncalled_file']).path, "7405745e675ea524d8e16e09b9bb749d")

    def test_portal_output_path_num_muts(self):
        self.assertEqualNumMutations([
            OFile.init_dict(self.res.output['portal_file']).path, 
            OFile.init_dict(self.res.output['uncalled_file']).path, 
            ], 
            OFile.init_dict(self.res.output['filtered_file']).path)

    def test_output_file_fields(self):
        self.assertMutFieldContains(
            OFile.init_dict(self.res.output['output_file']).path,
            "Tumor_Sample_Barcode", ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5"], containsAll = True)

    def test_portal_output_path_fields(self):
        self.assertMutFieldDoesntContain(
            OFile.init_dict(self.res.output['portal_file']).path,
            "Amino_Acid_Change", [""])

    def test_uncalled_output_path_fields(self):
        self.assertMutFieldDoesntContain(
            OFile.init_dict(self.res.output['uncalled_file']).path,
            "Amino_Acid_Change", [""])


class TestSamplesFilloutIndexBatch4Group(PlutoPreRunTestCase):
    
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
        ]

        sample_group2 = [
            {
                "sample_id": "Sample3",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "clinical",
                "prefilter": False,
                "maf_file": { "class": "File", "path": sample3_maf },
                "bam_file": { "class": "File", "path": sample3_bam }
            },
        ]

        sample_group3 = [
            {
                "sample_id": "Sample5",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "research",
                "prefilter": True,
                "maf_file": { "class": "File", "path": sample5_maf },
                "bam_file": { "class": "File", "path": sample5_bam }
            },
        ]

        sample_group4 = [
            {
                "sample_id": "Sample4",
                "normal_id": "FROZENPOOLEDNORMAL_IMPACT505_V2",
                "sample_type": "clinical",
                "prefilter": False,
                "maf_file": { "class": "File", "path": sample4_maf },
                "bam_file": { "class": "File", "path": sample4_bam }
            }
        ]

        self.input = {
            "sample_groups": [sample_group1, sample_group2, sample_group3, sample_group4],
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
        self.assertNumMutations(OFile.init_dict(self.res.output['output_file']).path, 222)

    def test_output_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['output_file']).path, "610d91b05a37b371a6b1b1615042dcc3")    

    def test_filtered_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['filtered_file']).path, 222)

    def test_filtered_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['filtered_file']).path, "ead4f5ba539e87e64eac2c5f3e02aba7")

    def test_portal_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['portal_file']).path, 113)

    def test_portal_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['portal_file']).path, "143f4d2e2eb0f09e9505e08c11baff03")

    def test_uncalled_file_num_muts(self):
        self.assertNumMutations(OFile.init_dict(self.res.output['uncalled_file']).path, 109)

    def test_uncalled_file_muts_hash(self):
        self.assertMutationsHash(OFile.init_dict(self.res.output['uncalled_file']).path, "788fe2049c3cb23d0232f046fcb98ba3")

    def test_portal_output_path_num_muts(self):
        self.assertEqualNumMutations([
            OFile.init_dict(self.res.output['portal_file']).path, 
            OFile.init_dict(self.res.output['uncalled_file']).path, 
            ], 
            OFile.init_dict(self.res.output['filtered_file']).path)

    def test_output_file_fields(self):
        self.assertMutFieldContains(
            OFile.init_dict(self.res.output['output_file']).path,
            "Tumor_Sample_Barcode", ["Sample1", "Sample2", "Sample3", "Sample4", "Sample5"], containsAll = True)

    def test_portal_output_path_fields(self):
        self.assertMutFieldDoesntContain(
            OFile.init_dict(self.res.output['portal_file']).path,
            "Amino_Acid_Change", [""])

    def test_uncalled_output_path_fields(self):
        self.assertMutFieldDoesntContain(
            OFile.init_dict(self.res.output['uncalled_file']).path,
            "Amino_Acid_Change", [""])



if __name__ == "__main__":
    unittest.main()
