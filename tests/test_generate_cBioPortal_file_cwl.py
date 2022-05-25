#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
unit tests for the generate_cbioPortal_file.cwl file
"""
import os
import sys
import unittest
from tempfile import TemporaryDirectory

THIS_DIR = os.path.dirname(os.path.realpath(__file__))
PARENT_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, PARENT_DIR)
from pluto.tools import run_cwl, write_table, CWLFile
from pluto.settings import DATA_SETS
sys.path.pop(0)

cwl_file = CWLFile('generate_cBioPortal_file.cwl')

class TestGenerateCbioFilesCWL(unittest.TestCase):
    def test_generate_meta_sample(self):
        """
        meta_clinical_sample.txt
        """
        input_json = {
        "subcommand": "meta_sample",
        "cancer_study_id": "cancer_study",
        "sample_data_filename": "data_clinical_sample.txt",
        "output_filename": "meta_clinical_sample.txt"
        }
        with TemporaryDirectory() as tmpdir:
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file
                )

            expected_output = {
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, "meta_clinical_sample.txt"),
                    'basename': "meta_clinical_sample.txt",
                    'class': 'File',
                    'checksum': 'sha1$14021d16e19aa53440f953aece0e66e41d09c7f5',
                    'size': 140,
                    'path': os.path.join(output_dir, "meta_clinical_sample.txt")
                    }
                }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

    def test_generate_data_clinical_patient(self):
        """
        data_clinical_patient.txt

        generate_cbioPortal_files.py \
        patient \
        --data-clinical-file "$(DATA_CLINICAL_FILE)" \
        --output "$(CBIO_CLINCIAL_PATIENT_DATA_FILE)"
        """
        data_clinical_file = os.path.join(DATA_SETS['Proj_08390_G']['INPUTS_DIR'], "Proj_08390_G_sample_data_clinical.txt")

        input_json = {
        "subcommand": "patient",
        "data_clinical_file": {
            "path": data_clinical_file,
            "class": "File"
            },
        "output_filename": "data_clinical_patient.txt"
        }
        with TemporaryDirectory() as tmpdir:
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file
                )

            expected_output = {
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'data_clinical_patient.txt'),
                    'basename': 'data_clinical_patient.txt',
                    'class': 'File',
                    'checksum': 'sha1$9417dcabddd6ab2cbe98167bccd9b9e4fa182562',
                    'size': 643,
                    'path': os.path.join(output_dir,'data_clinical_patient.txt')
                }
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

    def test_generate_data_clinical_sample(self):
        """
        # data_clinical_sample.txt

        generate_cbioPortal_files.py \
        sample \
        --data-clinical-file "$(DATA_CLINICAL_FILE)" \
        --sample-summary-file "$(SAMPLE_SUMMARY_FILE)" \
        --project-pi "$(PROJ_PI)" \
        --request-pi "$(REQUEST_PI)" \
        --output "$(CBIO_CLINICAL_SAMPLE_DATA_FILE)"
        """
        data_clinical_file = os.path.join(DATA_SETS['Proj_08390_G']['INPUTS_DIR'], "Proj_08390_G_sample_data_clinical.txt")
        sample_summary_file = os.path.join(DATA_SETS['Proj_08390_G']['QC_DIR'], "Proj_08390_G_SampleSummary.txt")

        input_json = {
        "subcommand": "sample",
        "data_clinical_file": {
            "path": data_clinical_file,
            "class": "File"
            },
        "sample_summary_file": {
            "path": sample_summary_file,
            "class": "File"
            },
        "output_filename": "data_clinical_sample.txt",
        "project_pi": "Dr. Jones",
        "request_pi": "Dr. Franklin"
        }

        with TemporaryDirectory() as tmpdir:
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file
                )

            expected_output = {
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'data_clinical_sample.txt'),
                    'basename': 'data_clinical_sample.txt',
                    'class': 'File',
                    'checksum': 'sha1$2a0c59593fa7726743b2fe46db9d955dbc625453',
                    'size': 7592,
                    'path': os.path.join(output_dir,'data_clinical_sample.txt')
                }
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

    def test_test_generate_data_clinical_sample_with_facets(self):
        """
        Test that the data clinical sample file is generated correctly when multiple Facets files are provided
        """
        data_clinical_file = os.path.join(DATA_SETS['Proj_08390_G']['INPUTS_DIR'], 'Proj_08390_G_sample_data_clinical.2.txt')
        facets_txt_file1 = os.path.join(DATA_SETS['Proj_08390_G']['FACETS_SUITE_DIR'], 'Sample46.txt')
        facets_txt_file2 = os.path.join(DATA_SETS['Proj_08390_G']['FACETS_SUITE_DIR'], 'Sample44.txt')
        input_json = {
        "subcommand": "sample",
        "data_clinical_file": {
            "path": data_clinical_file,
            "class": "File"
            },
        "output_filename": "data_clinical_sample.txt",
        "project_pi": "jonesd",
        "request_pi": "franklind",
        "facets_txt_files": [
            {
                "path": facets_txt_file1,
                "class": "File"
            },
            {
                "path": facets_txt_file2,
                "class": "File"
            }
        ]
        }
        with TemporaryDirectory() as tmpdir:
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file
                )

            expected_output = {
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'data_clinical_sample.txt'),
                    'basename': 'data_clinical_sample.txt',
                    'class': 'File',
                    'checksum': 'sha1$1a1aee93048facdfb3b25598e3b560a9b6b2856a',
                    'size': 1289,
                    'path': os.path.join(output_dir, 'data_clinical_sample.txt')
                }
            }
            self.assertDictEqual(output_json, expected_output)

            with open(os.path.join(output_dir, 'data_clinical_sample.txt')) as fin:
                lines = [ line.strip().split('\t') for line in fin ]
            expected_lines = [
            ['#SAMPLE_ID', 'IGO_ID', 'PATIENT_ID', 'COLLAB_ID', 'SAMPLE_TYPE', 'SAMPLE_CLASS', 'GENE_PANEL', 'ONCOTREE_CODE', 'SPECIMEN_PRESERVATION_TYPE', 'TISSUE_SITE', 'REQUEST_ID', 'PROJECT_ID', 'PIPELINE', 'PIPELINE_VERSION', 'PROJECT_PI', 'REQUEST_PI', 'genome_doubled', 'ASCN_PURITY', 'ASCN_PLOIDY', 'ASCN_VERSION', 'ASCN_WGD'],
            ['#SAMPLE_ID', 'IGO_ID', 'PATIENT_ID', 'COLLAB_ID', 'SAMPLE_TYPE', 'SAMPLE_CLASS', 'GENE_PANEL', 'ONCOTREE_CODE', 'SPECIMEN_PRESERVATION_TYPE', 'TISSUE_SITE', 'REQUEST_ID', 'PROJECT_ID', 'PIPELINE', 'PIPELINE_VERSION', 'PROJECT_PI', 'REQUEST_PI', 'genome_doubled', 'ASCN_PURITY', 'ASCN_PLOIDY', 'ASCN_VERSION', 'ASCN_WGD'],
            ['#STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'STRING', 'NUMBER', 'NUMBER', 'STRING', 'STRING'],
            ['#1', '1', '1', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '1', '1', '0', '1'],
            ['SAMPLE_ID', 'IGO_ID', 'PATIENT_ID', 'COLLAB_ID', 'SAMPLE_TYPE', 'SAMPLE_CLASS', 'GENE_PANEL', 'ONCOTREE_CODE', 'SPECIMEN_PRESERVATION_TYPE', 'TISSUE_SITE', 'REQUEST_ID', 'PROJECT_ID', 'PIPELINE', 'PIPELINE_VERSION', 'PROJECT_PI', 'REQUEST_PI', 'genome_doubled', 'ASCN_PURITY', 'ASCN_PLOIDY', 'ASCN_VERSION', 'ASCN_WGD'],
            ['Sample46', '08390_G_95', 'p_C_00001', 'COLLAB-01-T', 'Primary', 'Biopsy', 'IMPACT468+08390_Hg19', 'MEL', 'FFPE', '', '08390_G', '08390', 'roslin', '2.5.7', 'jonesd', 'franklind', 'FALSE', '0.36', '2.6', '0.5.14', 'no WGD'],
            ['Sample44', '08390_G_93', 'p_C_00002', 'COLLAB-01-T', 'Primary', 'Biopsy', 'IMPACT468+08390_Hg19', 'MEL', 'FFPE', '', '08390_G', '08390', 'roslin', '2.5.7', 'jonesd', 'franklind', 'FALSE', '0.51', '1.6', '0.5.14', 'no WGD']
            ]
            self.assertEqual(lines, expected_lines)



    def test_generate_meta_study(self):
        """
        # meta_study.txt
            generate_cbioPortal_files.py \
            study \
            --cancer-study-id "$(PROJ_ID)" \
            --name "$(PROJ_NAME)" \
            --short-name "$(PROJ_SHORT_NAME)" \
            --type-of-cancer "$(CANCER_TYPE)" \
            --description "$(PROJ_DESC)" \
            --output "$(CBIO_META_STUDY_FILE)" \
            $(EXTRA_GROUPS_STR)
        """
        input_json = {
        "subcommand": "study",
        "output_filename": "meta_study.txt",
        "cancer_study_id": "cancer_study",
        "name": "cancer_study",
        "short_name": "cancer_study",
        "type_of_cancer": "MEL",
        "description": "description",
        "extra_groups": "FOO1"
        }
        with TemporaryDirectory() as tmpdir:
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file
                )

            expected_output = {
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'meta_study.txt'),
                    'basename': 'meta_study.txt',
                    'class': 'File',
                    'checksum': 'sha1$9625b915f0eba999305026833fa8b32b6ebebaa0',
                    'size': 161,
                    'path': os.path.join(output_dir,'meta_study.txt')
                }
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

    def test_meta_clinical_patient(self):
        """
        # meta_clinical_patient.txt
        generate_cbioPortal_files.py \
        meta_patient \
        --cancer-study-id "$(PROJ_ID)" \
        --patient-data-filename "$(CBIO_CLINCIAL_PATIENT_DATA_FILENAME)" \
        --output "$(CBIO_CLINCAL_PATIENT_META_FILE)"
        """
        input_json = {
        "subcommand": "meta_patient",
        "output_filename": "meta_clinical_patient.txt",
        "cancer_study_id": "cancer_study",
        "patient_data_filename": "data_clinical_patient.txt"
        }
        with TemporaryDirectory() as tmpdir:
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file
                )

            expected_output = {
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'meta_clinical_patient.txt'),
                    'basename': 'meta_clinical_patient.txt',
                    'class': 'File',
                    'checksum': 'sha1$cae62ab4638ff2ff39b71a43b5bd996f8eea16ea',
                    'size': 142,
                    'path': os.path.join(output_dir,'meta_clinical_patient.txt')
                }
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

    def test_generate_meta_CNA(self):
        """
        # meta_CNA.txt
            generate_cbioPortal_files.py \
            meta_cna \
            --cancer-study-id "$(PROJ_ID)" \
            --cna-data-filename "$(CBIO_CNA_DATA_FILENAME)" \
            --output "$(CBIO_META_CNA_FILE)"
        """
        input_json = {
        "subcommand": "meta_cna",
        "output_filename": "meta_CNA.txt",
        "cancer_study_id": "cancer_study",
        "cna_data_filename": "data_CNA.txt"
        }
        with TemporaryDirectory() as tmpdir:
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file
                )

            expected_output = {
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'meta_CNA.txt'),
                    'basename': 'meta_CNA.txt',
                    'class': 'File',
                    'checksum': 'sha1$a0c50ba21af32710c6895201ec2ec74809f43fec',
                    'size': 270,
                    'path': os.path.join(output_dir,'meta_CNA.txt')
                }
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

    def test_generate_meta_fusion(self):
        """
        # meta_fusions.txt

        generate_cbioPortal_files.py \
        meta_fusion \
        --cancer-study-id "$(PROJ_ID)" \
        --fusion-data-filename "$(CBIO_FUSION_DATA_FILENAME)" \
        --output "$(CBIO_META_FUSIONS_FILE)"
        """
        input_json = {
        "subcommand": "meta_fusion",
        "output_filename": "meta_fusions.txt",
        "cancer_study_id": "cancer_study",
        "fusion_data_filename": "data_fusions.txt"
        }
        with TemporaryDirectory() as tmpdir:
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file
                )

            expected_output = {
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'meta_fusions.txt'),
                    'basename': 'meta_fusions.txt',
                    'class': 'File',
                    'checksum': 'sha1$5e71daac57615260e685b9f7184a86ddf0e3a6d4',
                    'size': 227,
                    'path': os.path.join(output_dir,'meta_fusions.txt')
                }
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

    def test_generate_meta_sv(self):
        """
        # meta_sv.txt

        generate_cbioPortal_files.py \
        meta_sv \
        --cancer-study-id "$(PROJ_ID)" \
        --sv-data-filename "$(CBIO_SV_DATA_FILENAME)" \
        --output "$(CBIO_META_SV_FILE)"
        """
        input_json = {
        "subcommand": "meta_sv",
        "output_filename": "meta_sv.txt",
        "cancer_study_id": "cancer_study",
        "sv_data_filename": "data_sv.txt"
        }
        with TemporaryDirectory() as tmpdir:
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file
                )

            expected_output = {
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'meta_sv.txt'),
                    'basename': 'meta_sv.txt',
                    'class': 'File',
                    'checksum': 'sha1$731d3681e1c7ff3207ed3a5475b599e52d0464e0',
                    'size': 257,
                    'path': os.path.join(output_dir,'meta_sv.txt')
                }
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

    def test_generate_meta_mutations_extended(self):
        """
        # meta_mutations_extended.txt

        generate_cbioPortal_files.py \
        meta_mutations \
        --cancer-study-id "$(PROJ_ID)" \
        --mutations-data-filename "$(CBIO_MUTATION_DATA_FILENAME)" \
        --output "$(CBIO_META_MUTATIONS_FILE)"
        """
        input_json = {
        "subcommand": "meta_mutations",
        "output_filename": "meta_mutations_extended.txt",
        "cancer_study_id": "cancer_study",
        "mutations_data_filename": "data_mutations_extended.txt"
        }
        with TemporaryDirectory() as tmpdir:
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file
                )

            expected_output = {
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'meta_mutations_extended.txt'),
                    'basename': 'meta_mutations_extended.txt',
                    'class': 'File',
                    'checksum': 'sha1$778c1cd5450ede4127023d65af177e0eb6b58db0',# 'sha1$d6681566b68ec2eba1c16369f6838ed52986b044',
                    'size': 270, #253,
                    'path': os.path.join(output_dir,'meta_mutations_extended.txt')
                }
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

    def test_generate_meta_segments(self):
        """
        # <project_id>_meta_cna_hg19_seg.txt

        generate_cbioPortal_files.py \
        meta_segments \
        --cancer-study-id "$(PROJ_ID)" \
        --output "$(CBIO_META_CNA_SEGMENTS_FILE)" \
        --segmented-data-file "$(CBIO_SEGMENT_DATA_FILENAME)"
        """
        input_json = {
        "subcommand": "meta_segments",
        "output_filename": "Proj_08390_G_meta_cna_hg19_seg.txt",
        "cancer_study_id": "cancer_study",
        "segmented_data_filename": "Proj_08390_G_data_cna_hg19.seg"
        }
        with TemporaryDirectory() as tmpdir:
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file
                )

            expected_output = {
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'Proj_08390_G_meta_cna_hg19_seg.txt'),
                    'basename': 'Proj_08390_G_meta_cna_hg19_seg.txt',
                    'class': 'File',
                    'checksum': 'sha1$72f05c56f8304f1e12f1d922ccfb89a3c8559660',
                    'size': 200,
                    'path': os.path.join(output_dir,'Proj_08390_G_meta_cna_hg19_seg.txt')
                }
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

    def test_generate_cases_all(self):
        """
        # cases_all.txt

        generate_cbioPortal_files.py \
        cases_all  \
        --cancer-study-id "$(PROJ_ID)" \
        --data-clinical-file "$(DATA_CLINICAL_FILE)" \
        --output "$(CBIO_CASES_ALL_FILE)"
        """
        data_clinical_file = os.path.join(DATA_SETS['Proj_08390_G']['INPUTS_DIR'], "Proj_08390_G_sample_data_clinical.txt")
        input_json = {
            "subcommand": "cases_all",
            "output_filename": "cases_all.txt",
            "data_clinical_file": {
                "path": data_clinical_file,
                "class": "File"
                },
            "cancer_study_id": "Proj_08390_G"
        }

        with TemporaryDirectory() as tmpdir:
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file
                )

            expected_output = {
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'cases_all.txt'),
                    'basename': 'cases_all.txt',
                    'class': 'File',
                    'checksum': 'sha1$b9e43289cec5603b0886b5e8507c8d019387c125',
                    'size': 616,
                    'path': os.path.join(output_dir,'cases_all.txt')
                }
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

    def test_generate_cases_cnaseq(self):
        """
        # cases_cnaseq.txt

        generate_cbioPortal_files.py \
        cases_cnaseq \
        --cancer-study-id "$(PROJ_ID)" \
        --data-clinical-file "$(DATA_CLINICAL_FILE)" \
        --output "$(CBIO_CASES_CNASEQ_FILE)"
        """
        data_clinical_file = os.path.join(DATA_SETS['Proj_08390_G']['INPUTS_DIR'], "Proj_08390_G_sample_data_clinical.txt")
        input_json = {
            "subcommand": "cases_cnaseq",
            "output_filename": "cases_cnaseq.txt",
            "data_clinical_file": {
                "path": data_clinical_file,
                "class": "File"
                },
            "cancer_study_id": "Proj_08390_G"
        }

        with TemporaryDirectory() as tmpdir:
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file
                )

            expected_output = {
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'cases_cnaseq.txt'),
                    'basename': 'cases_cnaseq.txt',
                    'class': 'File',
                    'checksum': 'sha1$b87e2da8dce0fddbadec348efe2986519b2a794b',
                    'size': 696,
                    'path': os.path.join(output_dir,'cases_cnaseq.txt')
                }
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

    def test_generate_cases_cna(self):
        """
        # cases_cna.txt

        generate_cbioPortal_files.py \
        cases_cna \
        --cancer-study-id "$(PROJ_ID)" \
        --data-clinical-file "$(DATA_CLINICAL_FILE)" \
        --output "$(CBIO_CASES_CNA_FILE)"
        """
        data_clinical_file = os.path.join(DATA_SETS['Proj_08390_G']['INPUTS_DIR'], "Proj_08390_G_sample_data_clinical.txt")
        input_json = {
            "subcommand": "cases_cna",
            "output_filename": "cases_cna.txt",
            "data_clinical_file": {
                "path": data_clinical_file,
                "class": "File"
                },
            "cancer_study_id": "Proj_08390_G"
        }

        with TemporaryDirectory() as tmpdir:
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file
                )

            expected_output = {
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'cases_cna.txt'),
                    'basename': 'cases_cna.txt',
                    'class': 'File',
                    'checksum': 'sha1$053481a8299e9430117f8e45e081aa7ec21033a6',
                    'size': 628,
                    'path': os.path.join(output_dir,'cases_cna.txt')
                }
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

    def test_generate_cases_sequenced(self):
        """
        # cases_sequenced.txt

        generate_cbioPortal_files.py \
        cases_sequenced \
        --cancer-study-id "$(PROJ_ID)" \
        --data-clinical-file "$(DATA_CLINICAL_FILE)" \
        --output "$(CBIO_CASES_SEQUENCED_FILE)"
        """
        data_clinical_file = os.path.join(DATA_SETS['Proj_08390_G']['INPUTS_DIR'], "Proj_08390_G_sample_data_clinical.txt")
        input_json = {
            "subcommand": "cases_sequenced",
            "output_filename": "cases_sequenced.txt",
            "data_clinical_file": {
                "path": data_clinical_file,
                "class": "File"
                },
            "cancer_study_id": "Proj_08390_G"
        }

        with TemporaryDirectory() as tmpdir:
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file
                )

            expected_output = {
                'output_file': {
                    'location': 'file://' + os.path.join(output_dir, 'cases_sequenced.txt'),
                    'basename': 'cases_sequenced.txt',
                    'class': 'File',
                    'checksum': 'sha1$ef9f5aef03c2527bf576470168660557ca1c7cc9',
                    'size': 641,
                    'path': os.path.join(output_dir,'cases_sequenced.txt')
                }
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

    def test_clean_facets_suite_cna_data(self):
        """
        Test that we can clean up the bad headers on some Facets Suite CNA files
        """
        with TemporaryDirectory() as tmpdir:
            cna_lines = [
            ['Hugo_Symbol', 'sample1', 'sample2_hisens'],
            ['ABL1', '3;1', '3;NA']
            ]
            input_file = write_table(tmpdir, "data_CNA.txt", cna_lines)

            input_json = {
                "subcommand": "clean_cna",
                "output_filename": "output.txt",
                "input_file": {
                    "path": input_file,
                    "class": "File"
                    }
            }

            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file
                )

            output_file = os.path.join(output_dir, 'output.txt')
            expected_output = {
                'output_file': {
                    'location': 'file://' + output_file,
                    'basename': 'output.txt',
                    'class': 'File',
                    'checksum': 'sha1$a0ce27a97eca74902ac63638e3d34b0ed77da554',
                    'size': 42,
                    'path': output_file
                }
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

            with open(output_file) as fin:
                lines = [ l for l in fin ]
            expected_lines = ['Hugo_Symbol\tsample1\tsample2\n', 'ABL1\t3;1\t3;NA\n']
            self.assertEqual(lines, expected_lines)

if __name__ == "__main__":
    unittest.main()
