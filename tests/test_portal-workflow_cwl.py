#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unit tests for the portal-workflow.cwl
"""
import os
import sys
import unittest
from tempfile import TemporaryDirectory

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import CWLFile, PlutoTestCase
from pluto.settings import DATA_SETS, KNOWN_FUSIONS_FILE, ENABLE_LARGE_TESTS
from pluto.serializer import OFile, ODir
sys.path.pop(0)


class TestPortalWorkflow(PlutoTestCase):
    cwl_file = CWLFile('portal-workflow.cwl')

    def test_run_worflow_one_maf(self):
        """
        Test that the workflow works correctly when run with a single maf
        """
        data_clinical_file = os.path.join(DATA_SETS['Proj_1']['INPUTS_DIR'], "Proj_1_sample_data_clinical.txt")
        sample_summary_file = os.path.join(DATA_SETS['Proj_1']['QC_DIR'], "Proj_1_SampleSummary.txt")
        self.input = {
            "project_id": "Proj_1",
            "project_name": "Proj_1",
            "project_short_name": "Proj_1",
            "project_description": "project",
            "project_pi": "Dr. Jones",
            "request_pi": "Dr. Franklin",
            "is_impact": True,
            "argos_version_string": "2.x",
            "cancer_type": "MEL",
            "cancer_study_identifier": 'Proj_1',
            "cbio_meta_cna_segments_filename": "Proj_1_meta_cna_hg19_seg.txt",
            "cbio_segment_data_filename": "Proj_1_data_cna_hg19.seg",
            "helix_filter_version": "20.06.1",
            "data_clinical_file": {
                "path": data_clinical_file,
                "class": "File"
            },
            "sample_summary_file": {
                "path": sample_summary_file,
                "class": "File"
            },
            "targets_list": {
                "path": DATA_SETS['Proj_1']["targets_list"],
                "class": "File"
            },
            "known_fusions_file": {
                "path": KNOWN_FUSIONS_FILE,
                "class": "File"
            },
            "mutation_maf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['MAF_DIR'], "Sample1.Sample2.muts.maf"),
                    "class": "File"
                }
            ],
            "mutation_svs_txt_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.portal.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_cncf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_DIR'], "Sample2.rg.md.abra.printreads__Sample1.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_seg_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_DIR'], "Sample2.rg.md.abra.printreads__Sample1.rg.md.abra.printreads_hisens.seg"),
                "class": "File"
                }
            ],
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'portal_case_list_dir': ODir(name='case_lists', items=[
                OFile(name='cases_all.txt', size=604, hash='b8abebd686a5a4e3897af9546ee3e9e5fee8cf25'),
                OFile(name='cases_cnaseq.txt', size=684, hash='3bc39070fa15b305c9b1f4a041b75c6a0ddbba9a'),
                OFile(name='cases_cna.txt', size=616, hash='8d3dc67168382830a270feec6a0efddf52cc5dd4'),
                OFile(name='cases_sequenced.txt', size=629, hash='188bde59b913c82b9a783d19c37b81c959b4696d')], dir=output_dir),
            'portal_clinical_patient_meta_file': OFile(
                name='meta_clinical_patient.txt', size=136, hash='bbfd617bded72d6e9f2071285ac5a7867b0ec6fb', dir=output_dir),
            'portal_cna_ascna_file': OFile(
                name='data_CNA.ascna.txt', size=6154, hash='e06f41b3ad538740519cec581bcbd5cd812e5f00', dir=output_dir),

            # NOTE: same as the merged_cna_file file in this case
            'portal_cna_data_file': OFile(
                name='data_CNA.txt', size=5355, hash='c1682f09406478fc404e40758b0c5e9c47dce5cb', dir=output_dir),

            'portal_data_clinical_patient_file': OFile(
                name='data_clinical_patient.txt', size=643, hash='9417dcabddd6ab2cbe98167bccd9b9e4fa182562', dir=output_dir),
            'portal_data_clinical_sample_file': OFile(
                name='data_clinical_sample.txt', size=7592, hash='2a0c59593fa7726743b2fe46db9d955dbc625453', dir=output_dir),
            'portal_fusions_data_file': OFile(
                name='data_fusions.txt', size=99, hash='c16f763b248813fcdde76f7486f1ddc4e9856038', dir=output_dir),
            'portal_hisens_segs': OFile(
                name='Proj_1_data_cna_hg19.seg', size=1322, hash='cae32a187d973441dd6e554e07ae81bebecb8980', dir=output_dir),
            'portal_meta_clinical_sample_file': OFile(
                name='meta_clinical_sample.txt', size=134, hash='29d7eda8ae439aaaa531b2d10fa5c03f943edf11', dir=output_dir),
            'portal_meta_cna_file': OFile(
                name='meta_CNA.txt', size=264, hash='1123609f24529c407b04b5dbd22efd6a453b3965', dir=output_dir),
            'portal_meta_cna_segments_file': OFile(
                name='Proj_1_meta_cna_hg19_seg.txt', size=188, hash='c100c7c6cfb7f67f991d356725abea6204e99d6b', dir=output_dir),
            'portal_meta_fusions_file': OFile(
                name='meta_fusions.txt', size=221, hash='5417138de92de1c35aa123c1e8800d710bb1f7cb', dir=output_dir),
            'portal_meta_mutations_extended_file': OFile(
                name='meta_mutations_extended.txt', size=264, hash='c1e0524b9ee612710b1921053bdb3f32120831ec', dir=output_dir),
            'portal_meta_study_file': OFile(
                name='meta_study.txt', size=134, hash='182c7c39315d7ce91cbb8d96f98134d676324cf6', dir=output_dir),
            'portal_muts_file': OFile(
                name='data_mutations_extended.txt', size=4766, hash='5c92b50bdd1d7c45af8fd180448a68926b8248e6', dir=output_dir),
            'portal_report': OFile(
                name='report.html', size=1016472, hash='4be1f9395bb83330dcffaecf76def4456db99a62', dir=output_dir)
        }

        self.maxDiff = None
        strip_related_keys = [
        ('basename', 'report.html', ['size', 'checksum'])
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)

        comments, mutations = self.load_mutations(os.path.join(output_dir,  'data_mutations_extended.txt'))
        self.assertEqual(len(mutations), 17)

        # load the data_CNA.txt file
        path = os.path.join(output_dir, 'data_CNA.txt') # renamed from the data_CNA.scna.txt file ...
        with open(path) as f:
            header = next(f)
        header_parts = header.split()
        expected_header_parts = ['Hugo_Symbol', 'Sample1']
        self.assertEqual(header_parts, expected_header_parts)

        path = os.path.join(output_dir, 'data_CNA.ascna.txt')
        with open(path) as f:
            header = next(f)
        header_parts = header.split()
        expected_header_parts = ['Hugo_Symbol', 'Sample1']
        self.assertEqual(header_parts, expected_header_parts)

        # Test that CNA file looks correct
        output_file = expected_output['portal_cna_data_file']['path']
        lines = self.read_table(output_file)
        self.assertEqual(lines[0], ['Hugo_Symbol', 'Sample1'])
        self.assertEqual(len(lines), 587)


    def test_run_worflow_two_mafs(self):
        """
        Test that the workflow works correctly when run with two maf files
        """
        data_clinical_file = os.path.join(DATA_SETS['Proj_1']['INPUTS_DIR'], "Proj_1_sample_data_clinical.txt")
        sample_summary_file = os.path.join(DATA_SETS['Proj_1']['QC_DIR'], "Proj_1_SampleSummary.txt")
        self.input = {
            "project_id": "Proj_1",
            "project_name": "Proj_1",
            "project_short_name": "Proj_1",
            "project_description": "project",
            "project_pi": "Dr. Jones",
            "request_pi": "Dr. Franklin",
            "is_impact": True,
            "argos_version_string": "2.x",
            "cancer_type": "MEL",
            "cancer_study_identifier": 'Proj_1',
            "cbio_meta_cna_segments_filename": "Proj_1_meta_cna_hg19_seg.txt",
            "cbio_segment_data_filename": "Proj_1_data_cna_hg19.seg",
            "helix_filter_version": "20.06.1",
            "data_clinical_file": {
                "path": data_clinical_file,
                "class": "File"
            },
            "sample_summary_file": {
                "path": sample_summary_file,
                "class": "File"
            },
            "targets_list": {
                "path": DATA_SETS['Proj_1']["targets_list"],
                "class": "File"
            },
            "known_fusions_file": {
                "path": KNOWN_FUSIONS_FILE,
                "class": "File"
            },
            "mutation_maf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['MAF_DIR'], "Sample1.Sample2.muts.maf"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['MAF_DIR'], "Sample4.Sample3.muts.maf"),
                    "class": "File"
                }
            ],
            "mutation_svs_txt_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.portal.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['MAF_DIR'], "Sample4.Sample3.svs.pass.vep.portal.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_cncf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_DIR'], "Sample2.rg.md.abra.printreads__Sample1.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_DIR'], "Sample3.rg.md.abra.printreads__Sample4.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_seg_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_DIR'], "Sample2.rg.md.abra.printreads__Sample1.rg.md.abra.printreads_hisens.seg"),
                "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_DIR'], "Sample3.rg.md.abra.printreads__Sample4.rg.md.abra.printreads_hisens.seg"),
                    "class": "File"
                }
            ]
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'portal_case_list_dir': ODir(name='case_lists', items=[
                OFile(name='cases_all.txt', size=604, hash='b8abebd686a5a4e3897af9546ee3e9e5fee8cf25'),
                OFile(name='cases_cnaseq.txt', size=684, hash='3bc39070fa15b305c9b1f4a041b75c6a0ddbba9a'),
                OFile(name='cases_cna.txt', size=616, hash='8d3dc67168382830a270feec6a0efddf52cc5dd4'),
                OFile(name='cases_sequenced.txt', size=629, hash='188bde59b913c82b9a783d19c37b81c959b4696d')], dir=output_dir),
            'portal_clinical_patient_meta_file': OFile(
                name='meta_clinical_patient.txt', size=136, hash='bbfd617bded72d6e9f2071285ac5a7867b0ec6fb', dir=output_dir),
            'portal_cna_ascna_file': OFile(
                name='data_CNA.ascna.txt', size=8769, hash='a6e509cd86a412ba398e3b756f468907dc92ae34', dir=output_dir),
            'portal_cna_data_file': OFile(
                name='data_CNA.txt', size=6764, hash='dd0e1181f412ddb30b4de4738f5b737e18a552be', dir=output_dir),
            'portal_data_clinical_patient_file': OFile(
                name='data_clinical_patient.txt', size=643, hash='9417dcabddd6ab2cbe98167bccd9b9e4fa182562', dir=output_dir),
            'portal_data_clinical_sample_file': OFile(
                name='data_clinical_sample.txt', size=7592, hash='2a0c59593fa7726743b2fe46db9d955dbc625453', dir=output_dir),
            'portal_fusions_data_file': OFile(
                name='data_fusions.txt', size=99, hash='c16f763b248813fcdde76f7486f1ddc4e9856038', dir=output_dir),
            'portal_hisens_segs': OFile(
                name='Proj_1_data_cna_hg19.seg', size=2571, hash='d93de81ab61012dcfa4cd9af0e2696b63ed26928', dir=output_dir),
            'portal_meta_clinical_sample_file': OFile(
                name='meta_clinical_sample.txt', size=134, hash='29d7eda8ae439aaaa531b2d10fa5c03f943edf11', dir=output_dir),
            'portal_meta_cna_file': OFile(
                name='meta_CNA.txt', size=264, hash='1123609f24529c407b04b5dbd22efd6a453b3965', dir=output_dir),
            'portal_meta_cna_segments_file': OFile(
                name='Proj_1_meta_cna_hg19_seg.txt', size=188, hash='c100c7c6cfb7f67f991d356725abea6204e99d6b', dir=output_dir),
            'portal_meta_fusions_file': OFile(
                name='meta_fusions.txt', size=221, hash='5417138de92de1c35aa123c1e8800d710bb1f7cb', dir=output_dir),
            'portal_meta_mutations_extended_file': OFile(
                name='meta_mutations_extended.txt', size=264, hash='c1e0524b9ee612710b1921053bdb3f32120831ec', dir=output_dir),
            'portal_meta_study_file': OFile(
                name='meta_study.txt', size=134, hash='182c7c39315d7ce91cbb8d96f98134d676324cf6', dir=output_dir),
            'portal_muts_file': OFile(
                name='data_mutations_extended.txt', size=6999, hash='5f225a2cc3c93c7a85e5a504dd648e00e0e9970c', dir=output_dir),
            'portal_report': OFile(
                name='report.html', size=1016698, hash='b71520700964c22eec0a0700bbe95018b2bf7bec', dir=output_dir)
            }

        strip_related_keys = [
        ('basename', 'report.html', ['size', 'checksum'])
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)

        comments, mutations = self.load_mutations(os.path.join(output_dir,  'data_mutations_extended.txt'))
        self.assertEqual(len(mutations), 27)

        # load the data_CNA.txt file
        path = os.path.join(output_dir, 'data_CNA.txt')
        with open(path) as f:
            header = next(f)
        header_parts = header.split()
        expected_header_parts = ['Hugo_Symbol', 'Sample1', 'Sample4']
        self.assertEqual(header_parts, expected_header_parts)

        path = os.path.join(output_dir, 'data_CNA.ascna.txt')
        with open(path) as f:
            header = next(f)
        header_parts = header.split()
        expected_header_parts = ['Hugo_Symbol', 'Sample1', 'Sample4']
        self.assertEqual(header_parts, expected_header_parts)

    def test_with_facets_txt(self):
        """
        Test that the workflow produces expected output when Facets Suite .txt files are added
        """
        # use reduced file for only these sample pairs
        # Sample45	Sample46
        # Sample43	Sample44
        data_clinical_file = os.path.join(DATA_SETS['Proj_1']['INPUTS_DIR'], "Proj_1_sample_data_clinical.2.txt")
        sample_summary_file = os.path.join(DATA_SETS['Proj_1']['QC_DIR'], "Proj_1_SampleSummary.txt")
        facets_txt_file1 = os.path.join(DATA_SETS['Proj_1']['FACETS_SUITE_DIR'], 'Sample46.txt')
        facets_txt_file1 = os.path.join(DATA_SETS['Proj_1']['FACETS_SUITE_DIR'], 'Sample44.txt')
        self.input = {
            "project_id": "Proj_1",
            "project_name": "Proj_1",
            "project_short_name": "Proj_1",
            "project_description": "project",
            "project_pi": "Dr. Jones",
            "request_pi": "Dr. Franklin",
            "is_impact": True,
            "argos_version_string": "2.x",
            "cancer_type": "MEL",
            "cancer_study_identifier": 'Proj_1',
            "cbio_meta_cna_segments_filename": "Proj_1_meta_cna_hg19_seg.txt",
            "cbio_segment_data_filename": "Proj_1_data_cna_hg19.seg",
            "helix_filter_version": "20.06.1",
            "data_clinical_file": {
                "path": data_clinical_file,
                "class": "File"
            },
            "sample_summary_file": {
                "path": sample_summary_file,
                "class": "File"
            },
            "targets_list": {
                "path": DATA_SETS['Proj_1']["targets_list"],
                "class": "File"
            },
            "known_fusions_file": {
                "path": KNOWN_FUSIONS_FILE,
                "class": "File"
            },
            "mutation_maf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['MAF_DIR'], "Sample46.Sample45.muts.maf"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['MAF_DIR'], "Sample44.Sample43.muts.maf"),
                    "class": "File"
                }
            ],
            "mutation_svs_txt_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['MAF_DIR'], "Sample46.Sample45.svs.pass.vep.portal.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['MAF_DIR'], "Sample44.Sample43.svs.pass.vep.portal.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_cncf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_DIR'], "Sample45.rg.md.abra.printreads__Sample46.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_DIR'], "Sample43.rg.md.abra.printreads__Sample44.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_seg_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_DIR'], "Sample45.rg.md.abra.printreads__Sample46.rg.md.abra.printreads_hisens.seg"),
                "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_DIR'], "Sample43.rg.md.abra.printreads__Sample44.rg.md.abra.printreads_hisens.seg"),
                "class": "File"
                }
            ],
            "facets_suite_txt_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_SUITE_DIR'], "Sample46.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_SUITE_DIR'], "Sample44.txt"),
                    "class": "File"
                },
            ],
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'portal_case_list_dir':
                ODir(name='case_lists', items=[
                OFile(name='cases_all.txt', size=194, hash='744b6aca3004f89bdb3f437ea87d4e3991394eaa'),
                OFile(name='cases_cnaseq.txt', size=274, hash='9f7b03a5d0a9b8911c5c76ed0e6fe7712ee2a21f'),
                OFile(name='cases_cna.txt', size=206, hash='cc34a33f2103c0a8f6630799fbc72508fd40aee4'),
                OFile(name='cases_sequenced.txt', size=219, hash='738b2c81d52fd202df1afda31d2750affd09b3f1')], dir=output_dir),
            'portal_clinical_patient_meta_file': OFile(
                name='meta_clinical_patient.txt', size=136, hash='bbfd617bded72d6e9f2071285ac5a7867b0ec6fb', dir=output_dir),
            'portal_cna_ascna_file': OFile(
                name='data_CNA.ascna.txt', size=8593, hash='b322943e957285327f68c7d6033032af07a47c65', dir=output_dir),
            'portal_cna_data_file': OFile(
                name='data_CNA.txt', size=6398, hash='b25e5d8ed7cf067448e96db712de542df8d564cd', dir=output_dir),
            'portal_data_clinical_patient_file': OFile(
                name='data_clinical_patient.txt', size=91, hash='e45f9904b4de2c75fd148798075af9f05848aa27', dir=output_dir),
            'portal_data_clinical_sample_file': OFile(
                name='data_clinical_sample.txt', size=1366, hash='b8d4e3c53b6bf9407a5d2b1f192fea3d51c020ae', dir=output_dir),
            'portal_fusions_data_file': OFile(
                name='data_fusions.txt', size=99, hash='c16f763b248813fcdde76f7486f1ddc4e9856038', dir=output_dir),
            'portal_hisens_segs': OFile(
                name='Proj_1_data_cna_hg19.seg', size=3321, hash='c591f3c35dd53df525c2e31b14942b2e5748e7d3', dir=output_dir),
            'portal_meta_clinical_sample_file': OFile(
                name='meta_clinical_sample.txt', size=134, hash='29d7eda8ae439aaaa531b2d10fa5c03f943edf11', dir=output_dir),
            'portal_meta_cna_file': OFile(
                name='meta_CNA.txt', size=264, hash='1123609f24529c407b04b5dbd22efd6a453b3965', dir=output_dir),
            'portal_meta_cna_segments_file': OFile(
                name='Proj_1_meta_cna_hg19_seg.txt', size=188, hash='c100c7c6cfb7f67f991d356725abea6204e99d6b', dir=output_dir),
            'portal_meta_fusions_file': OFile(
                name='meta_fusions.txt', size=221, hash='5417138de92de1c35aa123c1e8800d710bb1f7cb', dir=output_dir),
            'portal_meta_mutations_extended_file': OFile(
                name='meta_mutations_extended.txt', size=264, hash='c1e0524b9ee612710b1921053bdb3f32120831ec', dir=output_dir),
            'portal_meta_study_file': OFile(
                name='meta_study.txt', size=134, hash='182c7c39315d7ce91cbb8d96f98134d676324cf6', dir=output_dir),
            'portal_muts_file': OFile(
                name='data_mutations_extended.txt', size=4981, hash='ac0df0e945cab34a5dbd514d26e9125b4af4c52a', dir=output_dir),
            'portal_report': OFile(
                name='report.html', size=1012302, hash='f604e2f13c2e1ba63f5ec1611564cad103263a05', dir=output_dir)
            }

        self.maxDiff = None

        strip_related_keys = [
        ('basename', 'report.html', ['size', 'checksum'])
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)


        comments, mutations = self.load_mutations(os.path.join(output_dir,  'data_mutations_extended.txt'))
        self.assertEqual(len(mutations), 18)

        with open(os.path.join(output_dir,  'data_clinical_sample.txt')) as fin:
            lines = [ line.strip().split('\t') for line in fin ]

        expected_lines = [
        ['#SAMPLE_ID', 'IGO_ID', 'PATIENT_ID', 'COLLAB_ID', 'SAMPLE_TYPE', 'SAMPLE_CLASS', 'GENE_PANEL', 'ONCOTREE_CODE', 'SPECIMEN_PRESERVATION_TYPE', 'TISSUE_SITE', 'REQUEST_ID', 'PROJECT_ID', 'PIPELINE', 'PIPELINE_VERSION', 'SAMPLE_COVERAGE', 'PROJECT_PI', 'REQUEST_PI', 'genome_doubled', 'ASCN_PURITY', 'ASCN_PLOIDY', 'ASCN_VERSION', 'ASCN_WGD'],
        ['#SAMPLE_ID', 'IGO_ID', 'PATIENT_ID', 'COLLAB_ID', 'SAMPLE_TYPE', 'SAMPLE_CLASS', 'GENE_PANEL', 'ONCOTREE_CODE', 'SPECIMEN_PRESERVATION_TYPE', 'TISSUE_SITE', 'REQUEST_ID', 'PROJECT_ID', 'PIPELINE', 'PIPELINE_VERSION', 'SAMPLE_COVERAGE', 'PROJECT_PI', 'REQUEST_PI', 'genome_doubled', 'ASCN_PURITY', 'ASCN_PLOIDY', 'ASCN_VERSION', 'ASCN_WGD'],
        ['#STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'NUMBER', 'STRING', 'STRING', 'STRING', 'NUMBER', 'NUMBER', 'STRING', 'STRING'],
        ['#1', '1', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '1', '0', '1'],
        ['SAMPLE_ID', 'IGO_ID', 'PATIENT_ID', 'COLLAB_ID', 'SAMPLE_TYPE', 'SAMPLE_CLASS', 'GENE_PANEL', 'ONCOTREE_CODE', 'SPECIMEN_PRESERVATION_TYPE', 'TISSUE_SITE', 'REQUEST_ID', 'PROJECT_ID', 'PIPELINE', 'PIPELINE_VERSION', 'SAMPLE_COVERAGE', 'PROJECT_PI', 'REQUEST_PI', 'genome_doubled', 'ASCN_PURITY', 'ASCN_PLOIDY', 'ASCN_VERSION', 'ASCN_WGD'],
        ['Sample46', '08390_G_95', 'p_C_00001', 'COLLAB-01-T', 'Primary', 'Biopsy', 'IMPACT468+08390_Hg19', 'MEL', 'FFPE', '', '08390_G', '08390', 'roslin', '2.5.7', '108', 'Dr. Jones', 'Dr. Franklin', 'FALSE', '0.36', '2.6', '0.5.14', 'no WGD'],
        ['Sample44', '08390_G_93', 'p_C_00002', 'COLLAB-01-T', 'Primary', 'Biopsy', 'IMPACT468+08390_Hg19', 'MEL', 'FFPE', '', '08390_G', '08390', 'roslin', '2.5.7', '502', 'Dr. Jones', 'Dr. Franklin', 'FALSE', '0.51', '1.6', '0.5.14', 'no WGD']
        ]
        self.assertEqual(lines, expected_lines)

        # load the data_CNA.txt file
        path = os.path.join(output_dir, 'data_CNA.txt')
        with open(path) as f:
            header = next(f)
        header_parts = header.split()
        expected_header_parts = ['Hugo_Symbol', 'Sample44', 'Sample46']
        self.assertEqual(header_parts, expected_header_parts)

        path = os.path.join(output_dir, 'data_CNA.ascna.txt')
        with open(path) as f:
            header = next(f)
        header_parts = header.split()
        expected_header_parts = ['Hugo_Symbol', 'Sample44', 'Sample46']
        self.assertEqual(header_parts, expected_header_parts)

    def test_with_facets_txt_and_facets_mafs(self):
        """
        Test that the workflow produces expected output when Facets Suite .txt files are added and annotated .maf files from the Facets Suite workflow are used

        Also test that merged CNA files work correctly
        """
        lines1 = [
        ['Hugo_Symbol', 'Sample1', 'Sample2'],
        ["TAP1", "0", "0"],
        ["ERRFI1", "0", "0"],
        ["STK19", "", "0"],
        ]
        cna_file1 = self.write_table(self.tmpdir, filename = "cna1.txt", lines = lines1)
        # use reduced file for only these sample pairs
        # Sample45	Sample46
        # Sample43	Sample44
        data_clinical_file = os.path.join(DATA_SETS['Proj_1']['INPUTS_DIR'], "Proj_1_sample_data_clinical.2.txt")
        sample_summary_file = os.path.join(DATA_SETS['Proj_1']['QC_DIR'], "Proj_1_SampleSummary.txt")
        facets_txt_file1 = os.path.join(DATA_SETS['Proj_1']['FACETS_SUITE_DIR'], 'Sample46.txt')
        facets_txt_file1 = os.path.join(DATA_SETS['Proj_1']['FACETS_SUITE_DIR'], 'Sample44.txt')
        maf_file1 = os.path.join(DATA_SETS['Proj_1']['FACETS_SUITE_DIR'], "Sample46.Sample45_hisens.ccf.portal.maf")
        maf_file2 = os.path.join(DATA_SETS['Proj_1']['FACETS_SUITE_DIR'], "Sample44.Sample43_hisens.ccf.portal.maf")
        self.input = {
            "project_id": "Proj_1",
            "project_name": "Proj_1",
            "project_short_name": "Proj_1",
            "project_description": "project",
            "project_pi": "Dr. Jones",
            "request_pi": "Dr. Franklin",
            "is_impact": True,
            "argos_version_string": "2.x",
            "cancer_type": "MEL",
            "cancer_study_identifier": 'Proj_1',
            "cbio_meta_cna_segments_filename": "Proj_1_meta_cna_hg19_seg.txt",
            "cbio_segment_data_filename": "Proj_1_data_cna_hg19.seg",
            "helix_filter_version": "20.06.1",
            "data_clinical_file": {
                "path": data_clinical_file,
                "class": "File"
            },
            "sample_summary_file": {
                "path": sample_summary_file,
                "class": "File"
            },
            "targets_list": {
                "path": DATA_SETS['Proj_1']["targets_list"],
                "class": "File"
            },
            "known_fusions_file": {
                "path": KNOWN_FUSIONS_FILE,
                "class": "File"
            },
            "mutation_maf_files": [
                {
                    "path": maf_file1,
                    "class": "File"
                },
                {
                    "path": maf_file2,
                    "class": "File"
                }
            ],
            "mutation_svs_txt_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['MAF_DIR'], "Sample46.Sample45.svs.pass.vep.portal.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['MAF_DIR'], "Sample44.Sample43.svs.pass.vep.portal.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_cncf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_DIR'], "Sample45.rg.md.abra.printreads__Sample46.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_DIR'], "Sample43.rg.md.abra.printreads__Sample44.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_seg_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_DIR'], "Sample45.rg.md.abra.printreads__Sample46.rg.md.abra.printreads_hisens.seg"),
                "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_DIR'], "Sample43.rg.md.abra.printreads__Sample44.rg.md.abra.printreads_hisens.seg"),
                "class": "File"
                }
            ],
            "facets_suite_txt_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_SUITE_DIR'], "Sample46.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_SUITE_DIR'], "Sample44.txt"),
                    "class": "File"
                },
            ],
            "extra_cna_files": [{"class": "File", "path": cna_file1}]
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'portal_case_list_dir': ODir(name='case_lists', items=[
                OFile(name='cases_all.txt', size=194, hash='744b6aca3004f89bdb3f437ea87d4e3991394eaa'),
                OFile(name='cases_cnaseq.txt', size=274, hash='9f7b03a5d0a9b8911c5c76ed0e6fe7712ee2a21f'),
                OFile(name='cases_cna.txt', size=206, hash='cc34a33f2103c0a8f6630799fbc72508fd40aee4'),
                OFile(name='cases_sequenced.txt', size=219, hash='738b2c81d52fd202df1afda31d2750affd09b3f1')], dir=output_dir),
            'portal_clinical_patient_meta_file': OFile(
                name='meta_clinical_patient.txt', size=136, hash='bbfd617bded72d6e9f2071285ac5a7867b0ec6fb', dir=output_dir),
            'portal_cna_ascna_file': OFile(
                name='data_CNA.ascna.txt', size=8593, hash='b322943e957285327f68c7d6033032af07a47c65', dir=output_dir),
            'portal_cna_data_file': OFile(
                name='data_CNA.txt', size=9919, hash='ca32a20d522a1ab38816a800608a44b836ca96bf', dir=output_dir),
            'portal_data_clinical_patient_file': OFile(
                name='data_clinical_patient.txt', size=91, hash='e45f9904b4de2c75fd148798075af9f05848aa27', dir=output_dir),
            'portal_data_clinical_sample_file': OFile(
                name='data_clinical_sample.txt', size=1366, hash='b8d4e3c53b6bf9407a5d2b1f192fea3d51c020ae', dir=output_dir),
            'portal_fusions_data_file': OFile(
                name='data_fusions.txt', size=99, hash='c16f763b248813fcdde76f7486f1ddc4e9856038', dir=output_dir),
            'portal_hisens_segs': OFile(
                name='Proj_1_data_cna_hg19.seg', size=3321, hash='c591f3c35dd53df525c2e31b14942b2e5748e7d3', dir=output_dir),
            'portal_meta_clinical_sample_file': OFile(
                name='meta_clinical_sample.txt', size=134, hash='29d7eda8ae439aaaa531b2d10fa5c03f943edf11', dir=output_dir),
            'portal_meta_cna_file': OFile(
                name='meta_CNA.txt', size=264, hash='1123609f24529c407b04b5dbd22efd6a453b3965', dir=output_dir),
            'portal_meta_cna_segments_file': OFile(
                name='Proj_1_meta_cna_hg19_seg.txt', size=188, hash='c100c7c6cfb7f67f991d356725abea6204e99d6b', dir=output_dir),
            'portal_meta_fusions_file': OFile(
                name='meta_fusions.txt', size=221, hash='5417138de92de1c35aa123c1e8800d710bb1f7cb', dir=output_dir),
            'portal_meta_mutations_extended_file': OFile(
                name='meta_mutations_extended.txt', size=264, hash='c1e0524b9ee612710b1921053bdb3f32120831ec', dir=output_dir),
            'portal_meta_study_file': OFile(
                name='meta_study.txt', size=134, hash='182c7c39315d7ce91cbb8d96f98134d676324cf6', dir=output_dir),
            'portal_muts_file': OFile(
                name='data_mutations_extended.txt', size=5543, hash='5de6aee434d063fb0d544e4dcde68ca34bb4612b', dir=output_dir),
            'portal_report': OFile(
                name='report.html', size=1012302, hash='55da8083f3a82aac9f5fe0e90d8a7e7c9800ffe5', dir=output_dir)
            }

        self.maxDiff = None
        strip_related_keys = [
        ('basename', 'report.html', ['size', 'checksum'])
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)

        with open(os.path.join(output_dir,  'data_clinical_sample.txt')) as fin:
            lines = [ line.strip().split('\t') for line in fin ]

        expected_lines = [
        ['#SAMPLE_ID', 'IGO_ID', 'PATIENT_ID', 'COLLAB_ID', 'SAMPLE_TYPE', 'SAMPLE_CLASS', 'GENE_PANEL', 'ONCOTREE_CODE', 'SPECIMEN_PRESERVATION_TYPE', 'TISSUE_SITE', 'REQUEST_ID', 'PROJECT_ID', 'PIPELINE', 'PIPELINE_VERSION', 'SAMPLE_COVERAGE', 'PROJECT_PI', 'REQUEST_PI', 'genome_doubled', 'ASCN_PURITY', 'ASCN_PLOIDY', 'ASCN_VERSION', 'ASCN_WGD'],
        ['#SAMPLE_ID', 'IGO_ID', 'PATIENT_ID', 'COLLAB_ID', 'SAMPLE_TYPE', 'SAMPLE_CLASS', 'GENE_PANEL', 'ONCOTREE_CODE', 'SPECIMEN_PRESERVATION_TYPE', 'TISSUE_SITE', 'REQUEST_ID', 'PROJECT_ID', 'PIPELINE', 'PIPELINE_VERSION', 'SAMPLE_COVERAGE', 'PROJECT_PI', 'REQUEST_PI', 'genome_doubled', 'ASCN_PURITY', 'ASCN_PLOIDY', 'ASCN_VERSION', 'ASCN_WGD'],
        ['#STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'NUMBER', 'STRING', 'STRING', 'STRING', 'NUMBER', 'NUMBER', 'STRING', 'STRING'],
        ['#1', '1', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '1', '0', '1'],
        ['SAMPLE_ID', 'IGO_ID', 'PATIENT_ID', 'COLLAB_ID', 'SAMPLE_TYPE', 'SAMPLE_CLASS', 'GENE_PANEL', 'ONCOTREE_CODE', 'SPECIMEN_PRESERVATION_TYPE', 'TISSUE_SITE', 'REQUEST_ID', 'PROJECT_ID', 'PIPELINE', 'PIPELINE_VERSION', 'SAMPLE_COVERAGE', 'PROJECT_PI', 'REQUEST_PI', 'genome_doubled', 'ASCN_PURITY', 'ASCN_PLOIDY', 'ASCN_VERSION', 'ASCN_WGD'],
        ['Sample46', '08390_G_95', 'p_C_00001', 'COLLAB-01-T', 'Primary', 'Biopsy', 'IMPACT468+08390_Hg19', 'MEL', 'FFPE', '', '08390_G', '08390', 'roslin', '2.5.7', '108', 'Dr. Jones', 'Dr. Franklin', 'FALSE', '0.36', '2.6', '0.5.14', 'no WGD'],
        ['Sample44', '08390_G_93', 'p_C_00002', 'COLLAB-01-T', 'Primary', 'Biopsy', 'IMPACT468+08390_Hg19', 'MEL', 'FFPE', '', '08390_G', '08390', 'roslin', '2.5.7', '502', 'Dr. Jones', 'Dr. Franklin', 'FALSE', '0.51', '1.6', '0.5.14', 'no WGD']
        ]
        self.assertEqual(lines, expected_lines)

        comments, mutations = self.load_mutations(os.path.join(output_dir, 'data_mutations_extended.txt'))
        self.assertEqual(len(mutations), 18)

        colnames = mutations[0].keys()
        some_required_colnames = [
            "ASCN.TOTAL_COPY_NUMBER",
            "ASCN.MINOR_COPY_NUMBER",
            "ASCN.EXPECTED_ALT_COPIES",
            "ASCN.CCF_EXPECTED_COPIES",
            "ASCN.CCF_EXPECTED_COPIES_LOWER",
            "ASCN.CCF_EXPECTED_COPIES_UPPER",
            "ASCN.ASCN_METHOD",
            "ASCN.ASCN_INTEGER_COPY_NUMBER"
        ]
        for colname in some_required_colnames:
            self.assertTrue(colname in colnames)

        # load the data_CNA.txt file; this one includes the merged data !!
        path = os.path.join(output_dir, 'data_CNA.txt') # renamed from the data_CNA.scna.txt file ...
        with open(path) as f:
            header = next(f)
        header_parts = header.split()
        expected_header_parts = ['Hugo_Symbol', 'Sample44', 'Sample46', 'Sample1', 'Sample2']
        self.assertEqual(header_parts, expected_header_parts)

        path = os.path.join(output_dir, 'data_CNA.ascna.txt')
        with open(path) as f:
            header = next(f)
        header_parts = header.split()
        expected_header_parts = ['Hugo_Symbol', 'Sample44', 'Sample46']
        self.assertEqual(header_parts, expected_header_parts)

    @unittest.skipIf(ENABLE_LARGE_TESTS!=True, "is a large test")
    def test_with_mixed_mafs(self):
        """
        Test that the workflow produces expected output when both Argos maf files and Facets Suite maf files are used in the workflow
        """
        # use reduced file for only these sample pairs
        # Sample45	Sample46
        # Sample43	Sample44
        data_clinical_file = os.path.join(DATA_SETS['Proj_1']['INPUTS_DIR'], "Proj_1_sample_data_clinical.2.txt")
        maf_file1 = os.path.join(DATA_SETS['Proj_1']['MAF_DIR'], "Sample46.Sample45.muts.maf")
        maf_file2 = os.path.join(DATA_SETS['Proj_1']['FACETS_SUITE_DIR'], "Sample44.Sample43_hisens.ccf.portal.maf")
        self.input = {
            "project_id": "Proj_1",
            "project_name": "Proj_1",
            "project_short_name": "Proj_1",
            "project_description": "project",
            "project_pi": "Dr. Jones",
            "request_pi": "Dr. Franklin",
            "is_impact": True,
            "argos_version_string": "2.x",
            "cancer_type": "MEL",
            "cancer_study_identifier": 'Proj_1',
            "cbio_meta_cna_segments_filename": "Proj_1_meta_cna_hg19_seg.txt",
            "cbio_segment_data_filename": "Proj_1_data_cna_hg19.seg",
            "helix_filter_version": "20.06.1",
            "data_clinical_file": {
                "path": data_clinical_file,
                "class": "File"
            },
            "targets_list": {
                "path": DATA_SETS['Proj_1']["targets_list"],
                "class": "File"
            },
            "known_fusions_file": {
                "path": KNOWN_FUSIONS_FILE,
                "class": "File"
            },
            "mutation_maf_files": [
                {
                    "path": maf_file1,
                    "class": "File"
                },
                {
                    "path": maf_file2,
                    "class": "File"
                }
            ],
            "mutation_svs_txt_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['MAF_DIR'], "Sample46.Sample45.svs.pass.vep.portal.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['MAF_DIR'], "Sample44.Sample43.svs.pass.vep.portal.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_cncf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_DIR'], "Sample45.rg.md.abra.printreads__Sample46.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_DIR'], "Sample43.rg.md.abra.printreads__Sample44.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_seg_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_DIR'], "Sample45.rg.md.abra.printreads__Sample46.rg.md.abra.printreads_hisens.seg"),
                "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_1']['FACETS_DIR'], "Sample43.rg.md.abra.printreads__Sample44.rg.md.abra.printreads_hisens.seg"),
                "class": "File"
                }
            ]
        }

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'portal_case_list_dir': ODir(name='case_lists', items=[
                OFile(name='cases_all.txt', size=194, hash='744b6aca3004f89bdb3f437ea87d4e3991394eaa'),
                OFile(name='cases_cnaseq.txt', size=274, hash='9f7b03a5d0a9b8911c5c76ed0e6fe7712ee2a21f'),
                OFile(name='cases_cna.txt', size=206, hash='cc34a33f2103c0a8f6630799fbc72508fd40aee4'),
                OFile(name='cases_sequenced.txt', size=219, hash='738b2c81d52fd202df1afda31d2750affd09b3f1')], dir=output_dir),
            'portal_clinical_patient_meta_file': OFile(
                name='meta_clinical_patient.txt', size=136, hash='bbfd617bded72d6e9f2071285ac5a7867b0ec6fb', dir=output_dir),
            'portal_cna_ascna_file': OFile(
                name='data_CNA.ascna.txt', size=8593, hash='b322943e957285327f68c7d6033032af07a47c65', dir=output_dir),
            'portal_cna_data_file': OFile(
                name='data_CNA.txt', size=6398, hash='b25e5d8ed7cf067448e96db712de542df8d564cd', dir=output_dir),
            'portal_data_clinical_patient_file': OFile(
                name='data_clinical_patient.txt', size=91, hash='e45f9904b4de2c75fd148798075af9f05848aa27', dir=output_dir),
            'portal_data_clinical_sample_file': OFile(
                name='data_clinical_sample.txt', size=1015, hash='a53eac2ae9b8861019bc69d40fdef611db83498a', dir=output_dir),
            'portal_fusions_data_file': OFile(
                name='data_fusions.txt', size=99, hash='c16f763b248813fcdde76f7486f1ddc4e9856038', dir=output_dir),
            'portal_hisens_segs': OFile(
                name='Proj_1_data_cna_hg19.seg', size=3321, hash='c591f3c35dd53df525c2e31b14942b2e5748e7d3', dir=output_dir),
            'portal_meta_clinical_sample_file': OFile(
                name='meta_clinical_sample.txt', size=134, hash='29d7eda8ae439aaaa531b2d10fa5c03f943edf11', dir=output_dir),
            'portal_meta_cna_file': OFile(
                name='meta_CNA.txt', size=264, hash='1123609f24529c407b04b5dbd22efd6a453b3965', dir=output_dir),
            'portal_meta_cna_segments_file': OFile(
                name='Proj_1_meta_cna_hg19_seg.txt', size=188, hash='c100c7c6cfb7f67f991d356725abea6204e99d6b', dir=output_dir),
            'portal_meta_fusions_file': OFile(
                name='meta_fusions.txt', size=221, hash='5417138de92de1c35aa123c1e8800d710bb1f7cb', dir=output_dir),
            'portal_meta_mutations_extended_file': OFile(
                name='meta_mutations_extended.txt', size=264, hash='c1e0524b9ee612710b1921053bdb3f32120831ec', dir=output_dir),
            'portal_meta_study_file': OFile(
                name='meta_study.txt', size=134, hash='182c7c39315d7ce91cbb8d96f98134d676324cf6', dir=output_dir),
            'portal_muts_file': OFile(
                name='data_mutations_extended.txt', size=5673, hash='524e25fa067f00e856a7d48d4e9da5e90c4f157d', dir=output_dir),
            'portal_report': OFile(
                name='report.html', size=1012222, hash='e973fac6d2763fe3d2d7eafebaf0856638015d8b', dir=output_dir)
            }

        self.maxDiff = None
        strip_related_keys = [
        ('basename', 'report.html', ['size', 'checksum'])
        ]
        self.assertCWLDictEqual(output_json, expected_output, related_keys = strip_related_keys)

        comments, mutations = self.load_mutations(os.path.join(output_dir,  'data_mutations_extended.txt'))
        self.assertEqual(len(mutations), 18)

        colnames = mutations[0].keys()
        some_required_colnames = [
            "ASCN.TOTAL_COPY_NUMBER",
            "ASCN.MINOR_COPY_NUMBER",
            "ASCN.EXPECTED_ALT_COPIES",
            "ASCN.CCF_EXPECTED_COPIES",
            "ASCN.CCF_EXPECTED_COPIES_LOWER",
            "ASCN.CCF_EXPECTED_COPIES_UPPER",
            "ASCN.ASCN_METHOD",
            "ASCN.ASCN_INTEGER_COPY_NUMBER"
        ]

        for colname in some_required_colnames:
            self.assertTrue(colname in colnames, "Column label {} not present in the mutation file. Missing columns: {}".format(colname, [ c for c in some_required_colnames if c not in colnames ]))

        # load the data_CNA.txt file
        path = os.path.join(output_dir, 'data_CNA.txt')
        with open(path) as f:
            header = next(f)
        header_parts = header.split()
        expected_header_parts = ['Hugo_Symbol', 'Sample44', 'Sample46']
        self.assertEqual(header_parts, expected_header_parts)

        path = os.path.join(output_dir, 'data_CNA.ascna.txt')
        with open(path) as f:
            header = next(f)
        header_parts = header.split()
        expected_header_parts = ['Hugo_Symbol', 'Sample44', 'Sample46']
        self.assertEqual(header_parts, expected_header_parts)

if __name__ == "__main__":
    unittest.main()
