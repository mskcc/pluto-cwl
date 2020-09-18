#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unit tests for the portal-workflow.cwl
"""
import os
import json
import unittest
from tempfile import TemporaryDirectory, NamedTemporaryFile


# relative imports, from CLI and from parent project
if __name__ != "__main__":
    from .tools import run_command, load_mutations
    from .settings import CWL_DIR, CWL_ARGS, DATA_SETS, KNOWN_FUSIONS_FILE

if __name__ == "__main__":
    from tools import run_command, load_mutations
    from settings import CWL_DIR, CWL_ARGS, DATA_SETS, KNOWN_FUSIONS_FILE

cwl_file = os.path.join(CWL_DIR, 'portal-workflow.cwl')

class TestPortalWorkflow(unittest.TestCase):
    def test_run_worflow_one_maf(self):
        """
        Test that the workflow works correctly when run with a single maf
        """
        data_clinical_file = os.path.join(DATA_SETS['Proj_08390_G']['INPUTS_DIR'], "Proj_08390_G_sample_data_clinical.txt")
        sample_summary_file = os.path.join(DATA_SETS['Proj_08390_G']['QC_DIR'], "Proj_08390_G_SampleSummary.txt")
        input_json = {
            "project_id": "Proj_08390_G",
            "project_name": "Proj_08390_G",
            "project_short_name": "Proj_08390_G",
            "project_description": "project",
            "project_pi": "Dr. Jones",
            "request_pi": "Dr. Franklin",
            "is_impact": True,
            "argos_version_string": "2.x",
            "cancer_type": "MEL",
            "cancer_study_identifier": 'Proj_08390_G',
            "cbio_meta_cna_segments_filename": "Proj_08390_G_meta_cna_hg19_seg.txt",
            "cbio_segment_data_filename": "Proj_08390_G_data_cna_hg19.seg",
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
                "path": DATA_SETS['Proj_08390_G']["targets_list"],
                "class": "File"
            },
            "known_fusions_file": {
                "path": KNOWN_FUSIONS_FILE,
                "class": "File"
            },
            "mutation_maf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf"),
                    "class": "File"
                }
            ],
            "mutation_svs_txt_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.portal.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_cncf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample2.rg.md.abra.printreads__Sample1.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_seg_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample2.rg.md.abra.printreads__Sample1.rg.md.abra.printreads_hisens.seg"),
                "class": "File"
                }
            ],
        }

        with TemporaryDirectory() as tmpdir:

            input_json_file = os.path.join(tmpdir, "input.json")
            with open(input_json_file, "w") as json_out:
                json.dump(input_json, json_out)

            output_dir = os.path.join(tmpdir, "output")
            tmp_dir = os.path.join(tmpdir, "tmp")
            cache_dir = os.path.join(tmpdir, "cache")

            command = [
                "cwl-runner",
                *CWL_ARGS,
                "--outdir", output_dir,
                "--tmpdir-prefix", tmp_dir,
                "--cachedir", cache_dir,
                cwl_file, input_json_file
                ]

            returncode, proc_stdout, proc_stderr = run_command(command)

            if returncode != 0:
                print(proc_stderr)

            self.assertEqual(returncode, 0)

            output_json = json.loads(proc_stdout)

            expected_output = {
                    'portal_dir': {
                        'location': 'file://' + os.path.join(output_dir, 'portal'),
                        'path': os.path.join(output_dir, 'portal'),
                        'class': 'Directory',
                        'basename': 'portal',
                        'listing': [
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_clinical_sample.txt'),
                                'basename': 'meta_clinical_sample.txt',
                                'class': 'File',
                                'checksum': 'sha1$4c567d81c3b17a76c324fd3e2f73793a6e804f65',
                                'size': 140,
                                'path': os.path.join(output_dir, 'portal/meta_clinical_sample.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_clinical_patient.txt'),
                                'basename': 'data_clinical_patient.txt',
                                'class': 'File',
                                'checksum': 'sha1$9417dcabddd6ab2cbe98167bccd9b9e4fa182562',
                                'size': 643,
                                'path': os.path.join(output_dir, 'portal/data_clinical_patient.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_clinical_sample.txt'),
                                'basename': 'data_clinical_sample.txt',
                                'class': 'File',
                                'checksum': 'sha1$2a0c59593fa7726743b2fe46db9d955dbc625453',
                                'size': 7592,
                                'path': os.path.join(output_dir, 'portal/data_clinical_sample.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_study.txt'),
                                'basename': 'meta_study.txt',
                                'class': 'File',
                                'checksum': 'sha1$2b0a5fd1a97329adf7c3b1596c84cd6567059a95',
                                'size': 152,
                                'path': os.path.join(output_dir, 'portal/meta_study.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_clinical_patient.txt'),
                                'basename': 'meta_clinical_patient.txt',
                                'class': 'File',
                                'checksum': 'sha1$9cdc9a7e44a230c012f48b0236bdcf0bbc7de67f',
                                'size': 142,
                                'path': os.path.join(output_dir, 'portal/meta_clinical_patient.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_CNA.txt'),
                                'basename': 'meta_CNA.txt',
                                'class': 'File',
                                'checksum': 'sha1$a9bf16f6a0490b19e611e8814b85f7bf1d52417a',
                                'size': 270,
                                'path': os.path.join(output_dir, 'portal/meta_CNA.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_fusions.txt'),
                                'basename': 'meta_fusions.txt',
                                'class': 'File',
                                'checksum': 'sha1$77649e888bafc6a4ed61261d1c46d2f238e1c32b',
                                'size': 227,
                                'path': os.path.join(output_dir, 'portal/meta_fusions.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_mutations_extended.txt'),
                                'basename': 'meta_mutations_extended.txt',
                                'class': 'File',
                                'checksum': 'sha1$fd04fcd0129b35bb8b8aaef57b2efa16b8f42e1d',
                                'size': 253,
                                'path': os.path.join(output_dir, 'portal/meta_mutations_extended.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/Proj_08390_G_meta_cna_hg19_seg.txt'),
                                'basename': 'Proj_08390_G_meta_cna_hg19_seg.txt',
                                'class': 'File',
                                'checksum': 'sha1$59b54d3cd81acdd9fc21df1dc05a71cebfbfe11e',
                                'size': 200,
                                'path': os.path.join(output_dir, 'portal/Proj_08390_G_meta_cna_hg19_seg.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_CNA.txt'),
                                'basename': 'data_CNA.txt',
                                'class': 'File',
                                'checksum': 'sha1$931d82412733d7f93dd4117cd955f35e5dcbacc1',
                                'size': 5365,
                                'path': os.path.join(output_dir, 'portal/data_CNA.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_CNA.ascna.txt'),
                                'basename': 'data_CNA.ascna.txt',
                                'class': 'File',
                                'checksum': 'sha1$452d5ddef12a44693d5a98a05f5d300801734cfe',
                                'size': 6164,
                                'path': os.path.join(output_dir, 'portal/data_CNA.ascna.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_mutations_extended.txt'),
                                'basename': 'data_mutations_extended.txt',
                                'class': 'File',
                                'checksum': 'sha1$e713516cf04750a3e3f1ef932b1c7202d4b75bf2',
                                'size': 5106,
                                'path': os.path.join(output_dir, 'portal/data_mutations_extended.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/Proj_08390_G_data_cna_hg19.seg'),
                                'basename': 'Proj_08390_G_data_cna_hg19.seg',
                                'class': 'File',
                                'checksum': 'sha1$f0ebb82c34b6530447fa1e70b6dedcc039840d61',
                                'size': 1632,
                                'path': os.path.join(output_dir, 'portal/Proj_08390_G_data_cna_hg19.seg')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_fusions.txt'),
                                'basename': 'data_fusions.txt',
                                'class': 'File',
                                'checksum': 'sha1$c16f763b248813fcdde76f7486f1ddc4e9856038',
                                'size': 99,
                                'path': os.path.join(output_dir, 'portal/data_fusions.txt')
                            },
                            {
                                'class': 'Directory',
                                'basename': 'case_lists',
                                'location': 'file://' + os.path.join(output_dir, 'portal/case_lists'),
                                'path': os.path.join(output_dir, 'portal/case_lists'),
                                'listing': [
                                    {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_all.txt'),
                                    'basename': 'cases_all.txt',
                                    'class': 'File',
                                    'checksum': 'sha1$b9e43289cec5603b0886b5e8507c8d019387c125',
                                    'size': 616,
                                    'path': os.path.join(output_dir, 'portal/case_lists/cases_all.txt')},
                                    {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_cnaseq.txt'),
                                    'basename': 'cases_cnaseq.txt',
                                    'class': 'File',
                                    'checksum': 'sha1$b87e2da8dce0fddbadec348efe2986519b2a794b',
                                    'size': 696,
                                    'path': os.path.join(output_dir, 'portal/case_lists/cases_cnaseq.txt')},
                                    {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_cna.txt'),
                                    'basename': 'cases_cna.txt',
                                    'class': 'File',
                                    'checksum': 'sha1$053481a8299e9430117f8e45e081aa7ec21033a6',
                                    'size': 628,
                                    'path': os.path.join(output_dir, 'portal/case_lists/cases_cna.txt')},
                                    {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_sequenced.txt'),
                                    'basename': 'cases_sequenced.txt',
                                    'class': 'File',
                                    'checksum': 'sha1$ef9f5aef03c2527bf576470168660557ca1c7cc9',
                                    'size': 641,
                                    'path': os.path.join(output_dir, 'portal/case_lists/cases_sequenced.txt')}
                                ]
                            }
                        ]
                    }
                }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

            comments, mutations = load_mutations(os.path.join(output_dir, 'portal', 'data_mutations_extended.txt'))
            self.assertEqual(len(mutations), 17)

            # load the data_CNA.txt file
            path = os.path.join(output_dir, 'portal/data_CNA.txt') # renamed from the data_CNA.scna.txt file ...
            with open(path) as f:
                header = next(f)
            header_parts = header.split()
            expected_header_parts = ['Hugo_Symbol', 's_C_VJ7F47_P001_d']
            self.assertEqual(header_parts, expected_header_parts)

            path = os.path.join(output_dir, 'portal/data_CNA.ascna.txt')
            with open(path) as f:
                header = next(f)
            header_parts = header.split()
            expected_header_parts = ['Hugo_Symbol', 's_C_VJ7F47_P001_d']
            self.assertEqual(header_parts, expected_header_parts)

    def test_run_worflow_two_mafs(self):
        """
        Test that the workflow works correctly when run with two maf files
        """
        data_clinical_file = os.path.join(DATA_SETS['Proj_08390_G']['INPUTS_DIR'], "Proj_08390_G_sample_data_clinical.txt")
        sample_summary_file = os.path.join(DATA_SETS['Proj_08390_G']['QC_DIR'], "Proj_08390_G_SampleSummary.txt")
        input_json = {
            "project_id": "Proj_08390_G",
            "project_name": "Proj_08390_G",
            "project_short_name": "Proj_08390_G",
            "project_description": "project",
            "project_pi": "Dr. Jones",
            "request_pi": "Dr. Franklin",
            "is_impact": True,
            "argos_version_string": "2.x",
            "cancer_type": "MEL",
            "cancer_study_identifier": 'Proj_08390_G',
            "cbio_meta_cna_segments_filename": "Proj_08390_G_meta_cna_hg19_seg.txt",
            "cbio_segment_data_filename": "Proj_08390_G_data_cna_hg19.seg",
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
                "path": DATA_SETS['Proj_08390_G']["targets_list"],
                "class": "File"
            },
            "known_fusions_file": {
                "path": KNOWN_FUSIONS_FILE,
                "class": "File"
            },
            "mutation_maf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample4.Sample3.muts.maf"),
                    "class": "File"
                }
            ],
            "mutation_svs_txt_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.portal.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample4.Sample3.svs.pass.vep.portal.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_cncf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample2.rg.md.abra.printreads__Sample1.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample3.rg.md.abra.printreads__Sample4.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_seg_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample2.rg.md.abra.printreads__Sample1.rg.md.abra.printreads_hisens.seg"),
                "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample3.rg.md.abra.printreads__Sample4.rg.md.abra.printreads_hisens.seg"),
                    "class": "File"
                }
            ]
        }
        with TemporaryDirectory() as tmpdir:
            input_json_file = os.path.join(tmpdir, "input.json")
            with open(input_json_file, "w") as json_out:
                json.dump(input_json, json_out)

            output_dir = os.path.join(tmpdir, "output")
            tmp_dir = os.path.join(tmpdir, "tmp")
            cache_dir = os.path.join(tmpdir, "cache")

            command = [
                "cwl-runner",
                *CWL_ARGS,
                "--outdir", output_dir,
                "--tmpdir-prefix", tmp_dir,
                "--cachedir", cache_dir,
                cwl_file, input_json_file
                ]

            returncode, proc_stdout, proc_stderr = run_command(command)

            if returncode != 0:
                print(proc_stderr)

            self.assertEqual(returncode, 0)

            output_json = json.loads(proc_stdout)

            expected_output = {
                'portal_dir': {
                    'class': 'Directory',
                    'basename': 'portal',
                    'location': 'file://' + os.path.join(output_dir, 'portal'),
                    'path': os.path.join(output_dir, 'portal'),
                    'listing': [
                        {'location': 'file://' + os.path.join(output_dir, 'portal/meta_clinical_sample.txt'),
                        'basename': 'meta_clinical_sample.txt',
                        'class': 'File',
                        'checksum': 'sha1$4c567d81c3b17a76c324fd3e2f73793a6e804f65',
                        'size': 140,
                        'path': os.path.join(output_dir, 'portal/meta_clinical_sample.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/data_clinical_patient.txt'),
                        'basename': 'data_clinical_patient.txt',
                        'class': 'File',
                        'checksum': 'sha1$9417dcabddd6ab2cbe98167bccd9b9e4fa182562',
                        'size': 643,
                        'path': os.path.join(output_dir, 'portal/data_clinical_patient.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/data_clinical_sample.txt'),
                        'basename': 'data_clinical_sample.txt',
                        'class': 'File',
                        'checksum': 'sha1$2a0c59593fa7726743b2fe46db9d955dbc625453',
                        'size': 7592,
                        'path': os.path.join(output_dir, 'portal/data_clinical_sample.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/meta_study.txt'),
                        'basename': 'meta_study.txt',
                        'class': 'File',
                        'checksum': 'sha1$2b0a5fd1a97329adf7c3b1596c84cd6567059a95',
                        'size': 152,
                        'path': os.path.join(output_dir, 'portal/meta_study.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/meta_clinical_patient.txt'),
                        'basename': 'meta_clinical_patient.txt',
                        'class': 'File',
                        'checksum': 'sha1$9cdc9a7e44a230c012f48b0236bdcf0bbc7de67f',
                        'size': 142,
                        'path': os.path.join(output_dir, 'portal/meta_clinical_patient.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/meta_CNA.txt'),
                        'basename': 'meta_CNA.txt',
                        'class': 'File',
                        'checksum': 'sha1$a9bf16f6a0490b19e611e8814b85f7bf1d52417a',
                        'size': 270,
                        'path': os.path.join(output_dir, 'portal/meta_CNA.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/meta_fusions.txt'),
                        'basename': 'meta_fusions.txt',
                        'class': 'File',
                        'checksum': 'sha1$77649e888bafc6a4ed61261d1c46d2f238e1c32b',
                        'size': 227,
                        'path': os.path.join(output_dir, 'portal/meta_fusions.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/meta_mutations_extended.txt'),
                        'basename': 'meta_mutations_extended.txt',
                        'class': 'File',
                        'checksum': 'sha1$fd04fcd0129b35bb8b8aaef57b2efa16b8f42e1d',
                        'size': 253,
                        'path': os.path.join(output_dir, 'portal/meta_mutations_extended.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/Proj_08390_G_meta_cna_hg19_seg.txt'),
                        'basename': 'Proj_08390_G_meta_cna_hg19_seg.txt',
                        'class': 'File',
                        'checksum': 'sha1$59b54d3cd81acdd9fc21df1dc05a71cebfbfe11e',
                        'size': 200,
                        'path': os.path.join(output_dir, 'portal/Proj_08390_G_meta_cna_hg19_seg.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/data_CNA.txt'),
                        'basename': 'data_CNA.txt',
                        'class': 'File',
                        'checksum': 'sha1$09b4d944e50ea9d0e7567e04ce55b0f21d281255',
                        'size': 6784,
                        'path': os.path.join(output_dir, 'portal/data_CNA.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/data_CNA.ascna.txt'),
                        'basename': 'data_CNA.ascna.txt',
                        'class': 'File',
                        'checksum': 'sha1$d93ffe83137d9a77e2420b40ab3a2e0a1a5ad069',
                        'size': 8789,
                        'path': os.path.join(output_dir, 'portal/data_CNA.ascna.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/data_mutations_extended.txt'),
                        'basename': 'data_mutations_extended.txt',
                        'class': 'File',
                        'checksum': 'sha1$43469aa0f9125d3dca6217ee02641638c3a92e24',
                        'size': 7539,
                        'path': os.path.join(output_dir, 'portal/data_mutations_extended.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/Proj_08390_G_data_cna_hg19.seg'),
                        'basename': 'Proj_08390_G_data_cna_hg19.seg',
                        'class': 'File',
                        'checksum': 'sha1$f6a77b280c047a7e2082e3a09e8138f861790d3a',
                        'size': 3191,
                        'path': os.path.join(output_dir, 'portal/Proj_08390_G_data_cna_hg19.seg')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/data_fusions.txt'),
                        'basename': 'data_fusions.txt',
                        'class': 'File',
                        'checksum': 'sha1$c16f763b248813fcdde76f7486f1ddc4e9856038',
                        'size': 99,
                        'path': os.path.join(output_dir, 'portal/data_fusions.txt')},
                        {'class': 'Directory',
                        'basename': 'case_lists',
                        'listing': [{'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_all.txt'),
                        'basename': 'cases_all.txt',
                        'class': 'File',
                        'checksum': 'sha1$b9e43289cec5603b0886b5e8507c8d019387c125',
                        'size': 616,
                        'path': os.path.join(output_dir, 'portal/case_lists/cases_all.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_cnaseq.txt'),
                        'basename': 'cases_cnaseq.txt',
                        'class': 'File',
                        'checksum': 'sha1$b87e2da8dce0fddbadec348efe2986519b2a794b',
                        'size': 696,
                        'path': os.path.join(output_dir, 'portal/case_lists/cases_cnaseq.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_cna.txt'),
                        'basename': 'cases_cna.txt',
                        'class': 'File',
                        'checksum': 'sha1$053481a8299e9430117f8e45e081aa7ec21033a6',
                        'size': 628,
                        'path': os.path.join(output_dir, 'portal/case_lists/cases_cna.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_sequenced.txt'),
                        'basename': 'cases_sequenced.txt',
                        'class': 'File',
                        'checksum': 'sha1$ef9f5aef03c2527bf576470168660557ca1c7cc9',
                        'size': 641,
                        'path': os.path.join(output_dir, 'portal/case_lists/cases_sequenced.txt')}],
                        'location': 'file://' + os.path.join(output_dir, 'portal/case_lists'),
                        'path': os.path.join(output_dir, 'portal/case_lists')}
                    ]
                }
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

            comments, mutations = load_mutations(os.path.join(output_dir, 'portal', 'data_mutations_extended.txt'))
            self.assertEqual(len(mutations), 27)

            # load the data_CNA.txt file
            path = os.path.join(output_dir, 'portal/data_CNA.txt')
            with open(path) as f:
                header = next(f)
            header_parts = header.split()
            expected_header_parts = ['Hugo_Symbol', 's_C_VJ7F47_P001_d', 's_C_X50T9Y_P001_d']
            self.assertEqual(header_parts, expected_header_parts)

            path = os.path.join(output_dir, 'portal/data_CNA.ascna.txt')
            with open(path) as f:
                header = next(f)
            header_parts = header.split()
            expected_header_parts = ['Hugo_Symbol', 's_C_VJ7F47_P001_d', 's_C_X50T9Y_P001_d']
            self.assertEqual(header_parts, expected_header_parts)

    def test_with_facets_txt(self):
        """
        Test that the workflow produces expected output when Facets Suite .txt files are added
        """
        # use reduced file for only these sample pairs
        # Sample45	Sample46
        # Sample43	Sample44
        data_clinical_file = os.path.join(DATA_SETS['Proj_08390_G']['INPUTS_DIR'], "Proj_08390_G_sample_data_clinical.2.txt")
        sample_summary_file = os.path.join(DATA_SETS['Proj_08390_G']['QC_DIR'], "Proj_08390_G_SampleSummary.txt")
        facets_txt_file1 = os.path.join(DATA_SETS['Proj_08390_G']['FACETS_SUITE_DIR'], 'Sample46.txt')
        facets_txt_file1 = os.path.join(DATA_SETS['Proj_08390_G']['FACETS_SUITE_DIR'], 'Sample44.txt')
        input_json = {
            "project_id": "Proj_08390_G",
            "project_name": "Proj_08390_G",
            "project_short_name": "Proj_08390_G",
            "project_description": "project",
            "project_pi": "Dr. Jones",
            "request_pi": "Dr. Franklin",
            "is_impact": True,
            "argos_version_string": "2.x",
            "cancer_type": "MEL",
            "cancer_study_identifier": 'Proj_08390_G',
            "cbio_meta_cna_segments_filename": "Proj_08390_G_meta_cna_hg19_seg.txt",
            "cbio_segment_data_filename": "Proj_08390_G_data_cna_hg19.seg",
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
                "path": DATA_SETS['Proj_08390_G']["targets_list"],
                "class": "File"
            },
            "known_fusions_file": {
                "path": KNOWN_FUSIONS_FILE,
                "class": "File"
            },
            "mutation_maf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample46.Sample45.muts.maf"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample44.Sample43.muts.maf"),
                    "class": "File"
                }
            ],
            "mutation_svs_txt_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample46.Sample45.svs.pass.vep.portal.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample44.Sample43.svs.pass.vep.portal.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_cncf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample45.rg.md.abra.printreads__Sample46.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample43.rg.md.abra.printreads__Sample44.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_seg_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample45.rg.md.abra.printreads__Sample46.rg.md.abra.printreads_hisens.seg"),
                "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample43.rg.md.abra.printreads__Sample44.rg.md.abra.printreads_hisens.seg"),
                "class": "File"
                }
            ],
            "facets_suite_txt_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_SUITE_DIR'], "Sample46.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_SUITE_DIR'], "Sample44.txt"),
                    "class": "File"
                },
            ],
        }

        with TemporaryDirectory() as tmpdir:

            input_json_file = os.path.join(tmpdir, "input.json")
            with open(input_json_file, "w") as json_out:
                json.dump(input_json, json_out)

            output_dir = os.path.join(tmpdir, "output")
            tmp_dir = os.path.join(tmpdir, "tmp")
            cache_dir = os.path.join(tmpdir, "cache")

            command = [
                "cwl-runner",
                *CWL_ARGS,
                "--outdir", output_dir,
                "--tmpdir-prefix", tmp_dir,
                "--cachedir", cache_dir,
                cwl_file, input_json_file
                ]

            returncode, proc_stdout, proc_stderr = run_command(command)

            if returncode != 0:
                print(proc_stderr)

            self.assertEqual(returncode, 0)

            output_json = json.loads(proc_stdout)

            expected_output = {
                    'portal_dir': {
                        'location': 'file://' + os.path.join(output_dir, 'portal'),
                        'path': os.path.join(output_dir, 'portal'),
                        'class': 'Directory',
                        'basename': 'portal',
                        'listing': [
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_clinical_sample.txt'),
                                'basename': 'meta_clinical_sample.txt',
                                'class': 'File',
                                'checksum': 'sha1$4c567d81c3b17a76c324fd3e2f73793a6e804f65',
                                'size': 140,
                                'path': os.path.join(output_dir, 'portal/meta_clinical_sample.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_clinical_patient.txt'),
                                'basename': 'data_clinical_patient.txt',
                                'class': 'File',
                                'checksum': 'sha1$e45f9904b4de2c75fd148798075af9f05848aa27',
                                'size': 91,
                                'path': os.path.join(output_dir, 'portal/data_clinical_patient.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_clinical_sample.txt'),
                                'basename': 'data_clinical_sample.txt',
                                'class': 'File',
                                'checksum': 'sha1$b8d4e3c53b6bf9407a5d2b1f192fea3d51c020ae',
                                'size': 1366,
                                'path': os.path.join(output_dir, 'portal/data_clinical_sample.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_study.txt'),
                                'basename': 'meta_study.txt',
                                'class': 'File',
                                'checksum': 'sha1$2b0a5fd1a97329adf7c3b1596c84cd6567059a95',
                                'size': 152,
                                'path': os.path.join(output_dir, 'portal/meta_study.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_clinical_patient.txt'),
                                'basename': 'meta_clinical_patient.txt',
                                'class': 'File',
                                'checksum': 'sha1$9cdc9a7e44a230c012f48b0236bdcf0bbc7de67f',
                                'size': 142,
                                'path': os.path.join(output_dir, 'portal/meta_clinical_patient.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_CNA.txt'),
                                'basename': 'meta_CNA.txt',
                                'class': 'File',
                                'checksum': 'sha1$a9bf16f6a0490b19e611e8814b85f7bf1d52417a',
                                'size': 270,
                                'path': os.path.join(output_dir, 'portal/meta_CNA.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_fusions.txt'),
                                'basename': 'meta_fusions.txt',
                                'class': 'File',
                                'checksum': 'sha1$77649e888bafc6a4ed61261d1c46d2f238e1c32b',
                                'size': 227,
                                'path': os.path.join(output_dir, 'portal/meta_fusions.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_mutations_extended.txt'),
                                'basename': 'meta_mutations_extended.txt',
                                'class': 'File',
                                'checksum': 'sha1$fd04fcd0129b35bb8b8aaef57b2efa16b8f42e1d',
                                'size': 253,
                                'path': os.path.join(output_dir, 'portal/meta_mutations_extended.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/Proj_08390_G_meta_cna_hg19_seg.txt'),
                                'basename': 'Proj_08390_G_meta_cna_hg19_seg.txt',
                                'class': 'File',
                                'checksum': 'sha1$59b54d3cd81acdd9fc21df1dc05a71cebfbfe11e',
                                'size': 200,
                                'path': os.path.join(output_dir, 'portal/Proj_08390_G_meta_cna_hg19_seg.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_CNA.txt'),
                                'basename': 'data_CNA.txt',
                                'class': 'File',
                                'checksum': 'sha1$ec8541f7ecf6d843bf4d939a9b136ccc2f4578ef',
                                'size': 6416,
                                'path': os.path.join(output_dir, 'portal/data_CNA.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_CNA.ascna.txt'),
                                'basename': 'data_CNA.ascna.txt',
                                'class': 'File',
                                'checksum': 'sha1$dd7926b383b6f8dc9e5e7230b7b8eab8d16b99b4',
                                'size': 8611,
                                'path': os.path.join(output_dir, 'portal/data_CNA.ascna.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_mutations_extended.txt'),
                                'basename': 'data_mutations_extended.txt',
                                'class': 'File',
                                'checksum': 'sha1$4975ce31cdb3ca0e0895951f6331dff2998141c9',
                                'size': 5305,
                                'path': os.path.join(output_dir, 'portal/data_mutations_extended.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/Proj_08390_G_data_cna_hg19.seg'),
                                'basename': 'Proj_08390_G_data_cna_hg19.seg',
                                'class': 'File',
                                'checksum': 'sha1$ab238a36f1437145fe4ecf03283cfe6397a0ce27',
                                'size': 4032,
                                'path': os.path.join(output_dir, 'portal/Proj_08390_G_data_cna_hg19.seg')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_fusions.txt'),
                                'basename': 'data_fusions.txt',
                                'class': 'File',
                                'checksum': 'sha1$c16f763b248813fcdde76f7486f1ddc4e9856038',
                                'size': 99,
                                'path': os.path.join(output_dir, 'portal/data_fusions.txt')
                            },
                            {
                                'class': 'Directory',
                                'basename': 'case_lists',
                                'location': 'file://' + os.path.join(output_dir, 'portal/case_lists'),
                                'path': os.path.join(output_dir, 'portal/case_lists'),
                                'listing': [
                                    {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_all.txt'),
                                    'basename': 'cases_all.txt',
                                    'class': 'File',
                                    'checksum': 'sha1$11d9ddb461e7c26182e4a2fffeb1e974728c4b18',
                                    'size': 206,
                                    'path': os.path.join(output_dir, 'portal/case_lists/cases_all.txt')},
                                    {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_cnaseq.txt'),
                                    'basename': 'cases_cnaseq.txt',
                                    'class': 'File',
                                    'checksum': 'sha1$5deb87793fa5e1a2660bbb2a112edc48cade1a1e',
                                    'size': 286,
                                    'path': os.path.join(output_dir, 'portal/case_lists/cases_cnaseq.txt')},
                                    {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_cna.txt'),
                                    'basename': 'cases_cna.txt',
                                    'class': 'File',
                                    'checksum': 'sha1$fc8de96de4f948b49d5bb42d6023376038197873',
                                    'size': 218,
                                    'path': os.path.join(output_dir, 'portal/case_lists/cases_cna.txt')},
                                    {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_sequenced.txt'),
                                    'basename': 'cases_sequenced.txt',
                                    'class': 'File',
                                    'checksum': 'sha1$d1e36feb2024351123ad78f770718c12f66c3370',
                                    'size': 231,
                                    'path': os.path.join(output_dir, 'portal/case_lists/cases_sequenced.txt')}
                                ]
                            }
                        ]
                    }
                }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

            comments, mutations = load_mutations(os.path.join(output_dir, 'portal', 'data_mutations_extended.txt'))
            self.assertEqual(len(mutations), 18)

            with open(os.path.join(output_dir, 'portal', 'data_clinical_sample.txt')) as fin:
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
            path = os.path.join(output_dir, 'portal/data_CNA.txt')
            with open(path) as f:
                header = next(f)
            header_parts = header.split()
            expected_header_parts = ['Hugo_Symbol', 's_C_A11NF2_P001_d', 's_C_A5NEWD_P001_d']
            self.assertEqual(header_parts, expected_header_parts)

            path = os.path.join(output_dir, 'portal/data_CNA.ascna.txt')
            with open(path) as f:
                header = next(f)
            header_parts = header.split()
            expected_header_parts = ['Hugo_Symbol', 's_C_A11NF2_P001_d', 's_C_A5NEWD_P001_d']
            self.assertEqual(header_parts, expected_header_parts)

    def test_with_facets_txt_and_facets_mafs(self):
        """
        Test that the workflow produces expected output when Facets Suite .txt files are added and annotated .maf files from the Facets Suite workflow are used
        """
        # use reduced file for only these sample pairs
        # Sample45	Sample46
        # Sample43	Sample44
        data_clinical_file = os.path.join(DATA_SETS['Proj_08390_G']['INPUTS_DIR'], "Proj_08390_G_sample_data_clinical.2.txt")
        sample_summary_file = os.path.join(DATA_SETS['Proj_08390_G']['QC_DIR'], "Proj_08390_G_SampleSummary.txt")
        facets_txt_file1 = os.path.join(DATA_SETS['Proj_08390_G']['FACETS_SUITE_DIR'], 'Sample46.txt')
        facets_txt_file1 = os.path.join(DATA_SETS['Proj_08390_G']['FACETS_SUITE_DIR'], 'Sample44.txt')
        maf_file1 = os.path.join(DATA_SETS['Proj_08390_G']['FACETS_SUITE_DIR'], "Sample46.Sample45_hisens.ccf.portal.maf")
        maf_file2 = os.path.join(DATA_SETS['Proj_08390_G']['FACETS_SUITE_DIR'], "Sample44.Sample43_hisens.ccf.portal.maf")
        input_json = {
            "project_id": "Proj_08390_G",
            "project_name": "Proj_08390_G",
            "project_short_name": "Proj_08390_G",
            "project_description": "project",
            "project_pi": "Dr. Jones",
            "request_pi": "Dr. Franklin",
            "is_impact": True,
            "argos_version_string": "2.x",
            "cancer_type": "MEL",
            "cancer_study_identifier": 'Proj_08390_G',
            "cbio_meta_cna_segments_filename": "Proj_08390_G_meta_cna_hg19_seg.txt",
            "cbio_segment_data_filename": "Proj_08390_G_data_cna_hg19.seg",
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
                "path": DATA_SETS['Proj_08390_G']["targets_list"],
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
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample46.Sample45.svs.pass.vep.portal.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample44.Sample43.svs.pass.vep.portal.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_cncf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample45.rg.md.abra.printreads__Sample46.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample43.rg.md.abra.printreads__Sample44.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_seg_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample45.rg.md.abra.printreads__Sample46.rg.md.abra.printreads_hisens.seg"),
                "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample43.rg.md.abra.printreads__Sample44.rg.md.abra.printreads_hisens.seg"),
                "class": "File"
                }
            ],
            "facets_suite_txt_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_SUITE_DIR'], "Sample46.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_SUITE_DIR'], "Sample44.txt"),
                    "class": "File"
                },
            ],
        }

        with TemporaryDirectory() as tmpdir:

            input_json_file = os.path.join(tmpdir, "input.json")
            with open(input_json_file, "w") as json_out:
                json.dump(input_json, json_out)

            output_dir = os.path.join(tmpdir, "output")
            tmp_dir = os.path.join(tmpdir, "tmp")
            cache_dir = os.path.join(tmpdir, "cache")

            command = [
                "cwl-runner",
                *CWL_ARGS,
                "--debug",
                "--outdir", output_dir,
                "--tmpdir-prefix", tmp_dir,
                "--cachedir", cache_dir,
                cwl_file, input_json_file
                ]

            returncode, proc_stdout, proc_stderr = run_command(command)

            if returncode != 0:
                print(proc_stderr)

            self.assertEqual(returncode, 0)

            output_json = json.loads(proc_stdout)

            expected_output = {
                    'portal_dir': {
                        'location': 'file://' + os.path.join(output_dir, 'portal'),
                        'path': os.path.join(output_dir, 'portal'),
                        'class': 'Directory',
                        'basename': 'portal',
                        'listing': [
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_clinical_sample.txt'),
                                'basename': 'meta_clinical_sample.txt',
                                'class': 'File',
                                'checksum': 'sha1$4c567d81c3b17a76c324fd3e2f73793a6e804f65',
                                'size': 140,
                                'path': os.path.join(output_dir, 'portal/meta_clinical_sample.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_clinical_patient.txt'),
                                'basename': 'data_clinical_patient.txt',
                                'class': 'File',
                                'checksum': 'sha1$e45f9904b4de2c75fd148798075af9f05848aa27',
                                'size': 91,
                                'path': os.path.join(output_dir, 'portal/data_clinical_patient.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_clinical_sample.txt'),
                                'basename': 'data_clinical_sample.txt',
                                'class': 'File',
                                'checksum': 'sha1$b8d4e3c53b6bf9407a5d2b1f192fea3d51c020ae',
                                'size': 1366,
                                'path': os.path.join(output_dir, 'portal/data_clinical_sample.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_study.txt'),
                                'basename': 'meta_study.txt',
                                'class': 'File',
                                'checksum': 'sha1$2b0a5fd1a97329adf7c3b1596c84cd6567059a95',
                                'size': 152,
                                'path': os.path.join(output_dir, 'portal/meta_study.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_clinical_patient.txt'),
                                'basename': 'meta_clinical_patient.txt',
                                'class': 'File',
                                'checksum': 'sha1$9cdc9a7e44a230c012f48b0236bdcf0bbc7de67f',
                                'size': 142,
                                'path': os.path.join(output_dir, 'portal/meta_clinical_patient.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_CNA.txt'),
                                'basename': 'meta_CNA.txt',
                                'class': 'File',
                                'checksum': 'sha1$a9bf16f6a0490b19e611e8814b85f7bf1d52417a',
                                'size': 270,
                                'path': os.path.join(output_dir, 'portal/meta_CNA.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_fusions.txt'),
                                'basename': 'meta_fusions.txt',
                                'class': 'File',
                                'checksum': 'sha1$77649e888bafc6a4ed61261d1c46d2f238e1c32b',
                                'size': 227,
                                'path': os.path.join(output_dir, 'portal/meta_fusions.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_mutations_extended.txt'),
                                'basename': 'meta_mutations_extended.txt',
                                'class': 'File',
                                'checksum': 'sha1$fd04fcd0129b35bb8b8aaef57b2efa16b8f42e1d',
                                'size': 253,
                                'path': os.path.join(output_dir, 'portal/meta_mutations_extended.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/Proj_08390_G_meta_cna_hg19_seg.txt'),
                                'basename': 'Proj_08390_G_meta_cna_hg19_seg.txt',
                                'class': 'File',
                                'checksum': 'sha1$59b54d3cd81acdd9fc21df1dc05a71cebfbfe11e',
                                'size': 200,
                                'path': os.path.join(output_dir, 'portal/Proj_08390_G_meta_cna_hg19_seg.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_CNA.txt'),
                                'basename': 'data_CNA.txt',
                                'class': 'File',
                                'checksum': 'sha1$ec8541f7ecf6d843bf4d939a9b136ccc2f4578ef',
                                'size': 6416,
                                'path': os.path.join(output_dir, 'portal/data_CNA.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_CNA.ascna.txt'),
                                'basename': 'data_CNA.ascna.txt',
                                'class': 'File',
                                'checksum': 'sha1$dd7926b383b6f8dc9e5e7230b7b8eab8d16b99b4',
                                'size': 8611,
                                'path': os.path.join(output_dir, 'portal/data_CNA.ascna.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_mutations_extended.txt'),
                                'basename': 'data_mutations_extended.txt',
                                'class': 'File',
                                'checksum': 'sha1$5de6aee434d063fb0d544e4dcde68ca34bb4612b',
                                'size': 5543,
                                'path': os.path.join(output_dir, 'portal/data_mutations_extended.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/Proj_08390_G_data_cna_hg19.seg'),
                                'basename': 'Proj_08390_G_data_cna_hg19.seg',
                                'class': 'File',
                                'checksum': 'sha1$ab238a36f1437145fe4ecf03283cfe6397a0ce27',
                                'size': 4032,
                                'path': os.path.join(output_dir, 'portal/Proj_08390_G_data_cna_hg19.seg')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_fusions.txt'),
                                'basename': 'data_fusions.txt',
                                'class': 'File',
                                'checksum': 'sha1$c16f763b248813fcdde76f7486f1ddc4e9856038',
                                'size': 99,
                                'path': os.path.join(output_dir, 'portal/data_fusions.txt')
                            },
                            {
                                'class': 'Directory',
                                'basename': 'case_lists',
                                'location': 'file://' + os.path.join(output_dir, 'portal/case_lists'),
                                'path': os.path.join(output_dir, 'portal/case_lists'),
                                'listing': [
                                    {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_all.txt'),
                                    'basename': 'cases_all.txt',
                                    'class': 'File',
                                    'checksum': 'sha1$11d9ddb461e7c26182e4a2fffeb1e974728c4b18',
                                    'size': 206,
                                    'path': os.path.join(output_dir, 'portal/case_lists/cases_all.txt')},
                                    {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_cnaseq.txt'),
                                    'basename': 'cases_cnaseq.txt',
                                    'class': 'File',
                                    'checksum': 'sha1$5deb87793fa5e1a2660bbb2a112edc48cade1a1e',
                                    'size': 286,
                                    'path': os.path.join(output_dir, 'portal/case_lists/cases_cnaseq.txt')},
                                    {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_cna.txt'),
                                    'basename': 'cases_cna.txt',
                                    'class': 'File',
                                    'checksum': 'sha1$fc8de96de4f948b49d5bb42d6023376038197873',
                                    'size': 218,
                                    'path': os.path.join(output_dir, 'portal/case_lists/cases_cna.txt')},
                                    {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_sequenced.txt'),
                                    'basename': 'cases_sequenced.txt',
                                    'class': 'File',
                                    'checksum': 'sha1$d1e36feb2024351123ad78f770718c12f66c3370',
                                    'size': 231,
                                    'path': os.path.join(output_dir, 'portal/case_lists/cases_sequenced.txt')}
                                ]
                            }
                        ]
                    }
                }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

            with open(os.path.join(output_dir, 'portal', 'data_clinical_sample.txt')) as fin:
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

            comments, mutations = load_mutations(os.path.join(output_dir, 'portal', 'data_mutations_extended.txt'))
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

            # load the data_CNA.txt file
            path = os.path.join(output_dir, 'portal/data_CNA.txt') # renamed from the data_CNA.scna.txt file ...
            with open(path) as f:
                header = next(f)
            header_parts = header.split()
            expected_header_parts = ['Hugo_Symbol', 's_C_A11NF2_P001_d', 's_C_A5NEWD_P001_d']
            self.assertEqual(header_parts, expected_header_parts)

            path = os.path.join(output_dir, 'portal/data_CNA.ascna.txt')
            with open(path) as f:
                header = next(f)
            header_parts = header.split()
            expected_header_parts = ['Hugo_Symbol', 's_C_A11NF2_P001_d', 's_C_A5NEWD_P001_d']
            self.assertEqual(header_parts, expected_header_parts)

    def test_with_mixed_mafs(self):
        """
        Test that the workflow produces expected output when both Argos maf files and Facets Suite maf files are used in the workflow
        """
        # use reduced file for only these sample pairs
        # Sample45	Sample46
        # Sample43	Sample44
        data_clinical_file = os.path.join(DATA_SETS['Proj_08390_G']['INPUTS_DIR'], "Proj_08390_G_sample_data_clinical.2.txt")
        maf_file1 = os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample46.Sample45.muts.maf")
        maf_file2 = os.path.join(DATA_SETS['Proj_08390_G']['FACETS_SUITE_DIR'], "Sample44.Sample43_hisens.ccf.portal.maf")
        input_json = {
            "project_id": "Proj_08390_G",
            "project_name": "Proj_08390_G",
            "project_short_name": "Proj_08390_G",
            "project_description": "project",
            "project_pi": "Dr. Jones",
            "request_pi": "Dr. Franklin",
            "is_impact": True,
            "argos_version_string": "2.x",
            "cancer_type": "MEL",
            "cancer_study_identifier": 'Proj_08390_G',
            "cbio_meta_cna_segments_filename": "Proj_08390_G_meta_cna_hg19_seg.txt",
            "cbio_segment_data_filename": "Proj_08390_G_data_cna_hg19.seg",
            "helix_filter_version": "20.06.1",
            "data_clinical_file": {
                "path": data_clinical_file,
                "class": "File"
            },
            "targets_list": {
                "path": DATA_SETS['Proj_08390_G']["targets_list"],
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
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample46.Sample45.svs.pass.vep.portal.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample44.Sample43.svs.pass.vep.portal.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_cncf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample45.rg.md.abra.printreads__Sample46.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample43.rg.md.abra.printreads__Sample44.rg.md.abra.printreads_hisens.cncf.txt"),
                    "class": "File"
                }
            ],
            "facets_hisens_seg_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample45.rg.md.abra.printreads__Sample46.rg.md.abra.printreads_hisens.seg"),
                "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample43.rg.md.abra.printreads__Sample44.rg.md.abra.printreads_hisens.seg"),
                "class": "File"
                }
            ]
        }

        with TemporaryDirectory() as tmpdir:

            input_json_file = os.path.join(tmpdir, "input.json")
            with open(input_json_file, "w") as json_out:
                json.dump(input_json, json_out)

            output_dir = os.path.join(tmpdir, "output")
            tmp_dir = os.path.join(tmpdir, "tmp")
            cache_dir = os.path.join(tmpdir, "cache")

            command = [
                "cwl-runner",
                *CWL_ARGS,
                "--debug",
                "--outdir", output_dir,
                "--tmpdir-prefix", tmp_dir,
                "--cachedir", cache_dir,
                cwl_file, input_json_file
                ]

            returncode, proc_stdout, proc_stderr = run_command(command)

            if returncode != 0:
                print(proc_stderr)

            self.assertEqual(returncode, 0)

            output_json = json.loads(proc_stdout)

            expected_output = {
                    'portal_dir': {
                        'location': 'file://' + os.path.join(output_dir, 'portal'),
                        'path': os.path.join(output_dir, 'portal'),
                        'class': 'Directory',
                        'basename': 'portal',
                        'listing': [
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_clinical_sample.txt'),
                                'basename': 'meta_clinical_sample.txt',
                                'class': 'File',
                                'checksum': 'sha1$4c567d81c3b17a76c324fd3e2f73793a6e804f65',
                                'size': 140,
                                'path': os.path.join(output_dir, 'portal/meta_clinical_sample.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_clinical_patient.txt'),
                                'basename': 'data_clinical_patient.txt',
                                'class': 'File',
                                'checksum': 'sha1$e45f9904b4de2c75fd148798075af9f05848aa27',
                                'size': 91,
                                'path': os.path.join(output_dir, 'portal/data_clinical_patient.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_clinical_sample.txt'),
                                'basename': 'data_clinical_sample.txt',
                                'class': 'File',
                                'checksum': 'sha1$a53eac2ae9b8861019bc69d40fdef611db83498a',
                                'size': 1015,
                                'path': os.path.join(output_dir, 'portal/data_clinical_sample.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_study.txt'),
                                'basename': 'meta_study.txt',
                                'class': 'File',
                                'checksum': 'sha1$2b0a5fd1a97329adf7c3b1596c84cd6567059a95',
                                'size': 152,
                                'path': os.path.join(output_dir, 'portal/meta_study.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_clinical_patient.txt'),
                                'basename': 'meta_clinical_patient.txt',
                                'class': 'File',
                                'checksum': 'sha1$9cdc9a7e44a230c012f48b0236bdcf0bbc7de67f',
                                'size': 142,
                                'path': os.path.join(output_dir, 'portal/meta_clinical_patient.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_CNA.txt'),
                                'basename': 'meta_CNA.txt',
                                'class': 'File',
                                'checksum': 'sha1$a9bf16f6a0490b19e611e8814b85f7bf1d52417a',
                                'size': 270,
                                'path': os.path.join(output_dir, 'portal/meta_CNA.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_fusions.txt'),
                                'basename': 'meta_fusions.txt',
                                'class': 'File',
                                'checksum': 'sha1$77649e888bafc6a4ed61261d1c46d2f238e1c32b',
                                'size': 227,
                                'path': os.path.join(output_dir, 'portal/meta_fusions.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/meta_mutations_extended.txt'),
                                'basename': 'meta_mutations_extended.txt',
                                'class': 'File',
                                'checksum': 'sha1$fd04fcd0129b35bb8b8aaef57b2efa16b8f42e1d',
                                'size': 253,
                                'path': os.path.join(output_dir, 'portal/meta_mutations_extended.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/Proj_08390_G_meta_cna_hg19_seg.txt'),
                                'basename': 'Proj_08390_G_meta_cna_hg19_seg.txt',
                                'class': 'File',
                                'checksum': 'sha1$59b54d3cd81acdd9fc21df1dc05a71cebfbfe11e',
                                'size': 200,
                                'path': os.path.join(output_dir, 'portal/Proj_08390_G_meta_cna_hg19_seg.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_CNA.txt'),
                                'basename': 'data_CNA.txt',
                                'class': 'File',
                                'checksum': 'sha1$ec8541f7ecf6d843bf4d939a9b136ccc2f4578ef',
                                'size': 6416,
                                'path': os.path.join(output_dir, 'portal/data_CNA.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_CNA.ascna.txt'),
                                'basename': 'data_CNA.ascna.txt',
                                'class': 'File',
                                'checksum': 'sha1$dd7926b383b6f8dc9e5e7230b7b8eab8d16b99b4',
                                'size': 8611,
                                'path': os.path.join(output_dir, 'portal/data_CNA.ascna.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_mutations_extended.txt'),
                                'basename': 'data_mutations_extended.txt',
                                'class': 'File',
                                'checksum': 'sha1$4d45b75a3414a395aa6e2a3e561f487eda0b0e23',
                                'size': 5781,
                                'path': os.path.join(output_dir, 'portal/data_mutations_extended.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/Proj_08390_G_data_cna_hg19.seg'),
                                'basename': 'Proj_08390_G_data_cna_hg19.seg',
                                'class': 'File',
                                'checksum': 'sha1$ab238a36f1437145fe4ecf03283cfe6397a0ce27',
                                'size': 4032,
                                'path': os.path.join(output_dir, 'portal/Proj_08390_G_data_cna_hg19.seg')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_fusions.txt'),
                                'basename': 'data_fusions.txt',
                                'class': 'File',
                                'checksum': 'sha1$c16f763b248813fcdde76f7486f1ddc4e9856038',
                                'size': 99,
                                'path': os.path.join(output_dir, 'portal/data_fusions.txt')
                            },
                            {
                                'class': 'Directory',
                                'basename': 'case_lists',
                                'location': 'file://' + os.path.join(output_dir, 'portal/case_lists'),
                                'path': os.path.join(output_dir, 'portal/case_lists'),
                                'listing': [
                                    {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_all.txt'),
                                    'basename': 'cases_all.txt',
                                    'class': 'File',
                                    'checksum': 'sha1$11d9ddb461e7c26182e4a2fffeb1e974728c4b18',
                                    'size': 206,
                                    'path': os.path.join(output_dir, 'portal/case_lists/cases_all.txt')},
                                    {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_cnaseq.txt'),
                                    'basename': 'cases_cnaseq.txt',
                                    'class': 'File',
                                    'checksum': 'sha1$5deb87793fa5e1a2660bbb2a112edc48cade1a1e',
                                    'size': 286,
                                    'path': os.path.join(output_dir, 'portal/case_lists/cases_cnaseq.txt')},
                                    {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_cna.txt'),
                                    'basename': 'cases_cna.txt',
                                    'class': 'File',
                                    'checksum': 'sha1$fc8de96de4f948b49d5bb42d6023376038197873',
                                    'size': 218,
                                    'path': os.path.join(output_dir, 'portal/case_lists/cases_cna.txt')},
                                    {'location': 'file://' + os.path.join(output_dir, 'portal/case_lists/cases_sequenced.txt'),
                                    'basename': 'cases_sequenced.txt',
                                    'class': 'File',
                                    'checksum': 'sha1$d1e36feb2024351123ad78f770718c12f66c3370',
                                    'size': 231,
                                    'path': os.path.join(output_dir, 'portal/case_lists/cases_sequenced.txt')}
                                ]
                            }
                        ]
                    }
                }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

            comments, mutations = load_mutations(os.path.join(output_dir, 'portal', 'data_mutations_extended.txt'))
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
            path = os.path.join(output_dir, 'portal/data_CNA.txt')
            with open(path) as f:
                header = next(f)
            header_parts = header.split()
            expected_header_parts = ['Hugo_Symbol', 's_C_A11NF2_P001_d', 's_C_A5NEWD_P001_d']
            self.assertEqual(header_parts, expected_header_parts)

            path = os.path.join(output_dir, 'portal/data_CNA.ascna.txt')
            with open(path) as f:
                header = next(f)
            header_parts = header.split()
            expected_header_parts = ['Hugo_Symbol', 's_C_A11NF2_P001_d', 's_C_A5NEWD_P001_d']
            self.assertEqual(header_parts, expected_header_parts)

if __name__ == "__main__":
    unittest.main()
