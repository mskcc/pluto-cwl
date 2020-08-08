#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unit tests for the workflow_with_facets.cwl
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

cwl_file = os.path.join(CWL_DIR, 'workflow_with_facets.cwl')

class TestWorkflow(unittest.TestCase):
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
            "is_impact": "True",
            "argos_version_string": "2.x",
            "cancer_type": "MEL",
            "cancer_study_identifier": 'Proj_08390_G',
            "analysis_gene_cna_filename": "Proj_08390_G.gene.cna.txt",
            "analysis_mutations_filename": "Proj_08390_G.muts.maf",
            "analysis_segment_cna_filename": "Proj_08390_G.seg.cna.txt",
            "analysis_sv_filename": "Proj_08390_G.svs.maf",
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
            "mutation_svs_txt_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.portal.txt"),
                    "class": "File"
                }
            ],
            "mutation_svs_maf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.maf"),
                    "class": "File"
                }
            ],
            "pairs": [
                {
                    "pair_maf": {
                        "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf"),
                        "class": "File"
                    },
                    "snp_pileup": {
                         "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample2.rg.md.abra.printreads__Sample1.rg.md.abra.printreads.dat.gz"),
                         "class": "File"
                    },
                    "pair_id": "Sample1.Sample2",
                    "tumor_id": "Sample1",
                    "normal_id": "Sample2"
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
                'analysis_dir': {
                    'class': 'Directory',
                    'basename': 'analysis',
                    'location': 'file://' + os.path.join(output_dir, 'analysis'),
                    'path': os.path.join(output_dir, 'analysis'),
                    'listing': [
                        {
                            'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.gene.cna.txt'),
                            'basename': 'Proj_08390_G.gene.cna.txt',
                            'class': 'File',
                            'checksum': 'sha1$890129aa22c2e6b0e6568d4799c628748fdc5139',
                            'size': 85269,
                            'path': os.path.join(output_dir, 'analysis/Proj_08390_G.gene.cna.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf'),
                            'basename': 'Proj_08390_G.muts.maf',
                            'class': 'File',
                            'checksum': 'sha1$ab415a160a0827449e5f084568058b1a0893c18b',
                            'size': 38093,
                            'path': os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.seg.cna.txt'),
                            'basename': 'Proj_08390_G.seg.cna.txt',
                            'class': 'File',
                            'checksum': 'sha1$76e2c8a0a9c4200e26b8fc4a60701bd7b1b86560',
                            'size': 1459,
                            'path': os.path.join(output_dir, 'analysis/Proj_08390_G.seg.cna.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.svs.maf'),
                            'basename': 'Proj_08390_G.svs.maf',
                            'class': 'File',
                            'checksum': 'sha1$df420706bb5b772a79317843c0a01a3c88a9571d',
                            'size': 23603,
                            'path': os.path.join(output_dir, 'analysis/Proj_08390_G.svs.maf')
                            }
                        ]
                    },
                    'facets_dir': {
                        'location': 'file://' + os.path.join(output_dir, 'facets'),
                        'basename': 'facets',
                        'class': 'Directory',
                        'path': os.path.join(output_dir, 'facets'),
                        'listing': [{'basename': 'Sample1.Sample2',
                              'class': 'Directory',
                              'listing': [{'basename': 'Sample1.Sample2_hisens.ccf.portal.maf',
                                           'checksum': 'sha1$fa311e483303fbda7778622af7be57c7d62aac76',
                                           'class': 'File',
                                           'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.Sample2_hisens.ccf.portal.maf'),
                                           'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.Sample2_hisens.ccf.portal.maf'),
                                           'size': 19822273},
                                          {'basename': 'Sample1.arm_level.txt',
                                           'checksum': 'sha1$4f1d6620f3c613e40a147338739d3267e891dda7',
                                           'class': 'File',
                                           'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.arm_level.txt'),
                                           'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.arm_level.txt'),
                                           'size': 1986},
                                          {'basename': 'Sample1.txt',
                                           'checksum': 'sha1$8abe9bed5d1cb807b451c0647f658394d3f96708',
                                           'class': 'File',
                                           'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.txt'),
                                           'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.txt'),
                                           'size': 560},
                                          {'basename': 'Sample1.gene_level.txt',
                                           'checksum': 'sha1$6b31bff7758586ec39c94a0c2ddb8201eda5e0ce',
                                           'class': 'File',
                                           'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.gene_level.txt'),
                                           'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.gene_level.txt'),
                                           'size': 252793},
                                          {'basename': 'Sample1_hisens.cncf.txt',
                                           'checksum': 'sha1$a183b43abb9b7bb8f67ec6a0411539643d1275fe',
                                           'class': 'File',
                                           'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.cncf.txt'),
                                           'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.cncf.txt'),
                                           'size': 5774},
                                          {'basename': 'Sample1_hisens.rds',
                                           'checksum': 'sha1$d46031ce35efd7341146cac0db809c8cb75c97c2',
                                           'class': 'File',
                                           'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.rds'),
                                           'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.rds'),
                                           'size': 550930},
                                          {'basename': 'Sample1_hisens.seg',
                                           'checksum': 'sha1$1e0f1ad1cc751d53f02ad19717921027c3c96dfb',
                                           'class': 'File',
                                           'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.seg'),
                                           'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.seg'),
                                           'size': 1786},
                                          {'basename': 'Sample1_purity.rds',
                                           'checksum': 'sha1$265c7b94df391ae12e3bc1cc3ec4e763629a4f0a',
                                           'class': 'File',
                                           'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_purity.rds'),
                                           'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_purity.rds'),
                                           'size': 550388},
                                          {'basename': 'Sample1_purity.seg',
                                           'checksum': 'sha1$36904773ac5f194eac7f816a5fa9b055931bb48c',
                                           'class': 'File',
                                           'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_purity.seg'),
                                           'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_purity.seg'),
                                           'size': 1284},
                                          {'basename': 'Sample1.qc.txt',
                                           'checksum': 'sha1$6fdefe96a23c0bc52050c46a16625fca603e4eb6',
                                           'class': 'File',
                                           'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.qc.txt'),
                                           'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.qc.txt'),
                                           'size': 1340}],
                              'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2'),
                              'path': os.path.join(output_dir,'facets/Sample1.Sample2')
                                }]
                    },
                    'facets_failed_pairs': [],
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
                                'checksum': 'sha1$c285f8df3a2103e651023804ac2189ccc85dd3ef',
                                'size': 8550,
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
                                'checksum': 'sha1$e841a147242025dbb5e4980e3273be90469a2c34',
                                'size': 5326,
                                'path': os.path.join(output_dir, 'portal/data_CNA.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_CNA.ascna.txt'),
                                'basename': 'data_CNA.ascna.txt',
                                'class': 'File',
                                'checksum': 'sha1$9ee4be1ca109c92c62d0b7a6f5e976bf053646c9',
                                'size': 6174,
                                'path': os.path.join(output_dir, 'portal/data_CNA.ascna.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/data_mutations_extended.txt'),
                                'basename': 'data_mutations_extended.txt',
                                'class': 'File',
                                'checksum': 'sha1$ae78549d387625aa58358f3e52928d977273e4b4',
                                'size': 5629,
                                'path': os.path.join(output_dir, 'portal/data_mutations_extended.txt')
                            },
                            {
                                'location': 'file://' + os.path.join(output_dir, 'portal/Proj_08390_G_data_cna_hg19.seg'),
                                'basename': 'Proj_08390_G_data_cna_hg19.seg',
                                'class': 'File',
                                'checksum': 'sha1$76e2c8a0a9c4200e26b8fc4a60701bd7b1b86560',
                                'size': 1459,
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
            comments, mutations = load_mutations(os.path.join(output_dir, 'analysis', 'Proj_08390_G.muts.maf'))
            self.assertEqual(len(mutations), 22)
            comments, mutations = load_mutations(os.path.join(output_dir, 'portal', 'data_mutations_extended.txt'))
            self.assertEqual(len(mutations), 17)

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
            "is_impact": "True",
            "argos_version_string": "2.x",
            "cancer_type": "MEL",
            "cancer_study_identifier": 'Proj_08390_G',
            "analysis_gene_cna_filename": "Proj_08390_G.gene.cna.txt",
            "analysis_mutations_filename": "Proj_08390_G.muts.maf",
            "analysis_segment_cna_filename": "Proj_08390_G.seg.cna.txt",
            "analysis_sv_filename": "Proj_08390_G.svs.maf",
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
            "mutation_svs_maf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.maf"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample4.Sample3.svs.pass.vep.maf"),
                    "class": "File"
                }
            ],
            "pairs": [
                {
                    "pair_maf": {
                        "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf"),
                        "class": "File"
                    },
                    "snp_pileup": {
                         "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample2.rg.md.abra.printreads__Sample1.rg.md.abra.printreads.dat.gz"),
                         "class": "File"
                    },
                    "pair_id": "Sample1.Sample2",
                    "tumor_id": "Sample1",
                    "normal_id": "Sample2"
                },
                {
                    "pair_maf": {
                        "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample4.Sample3.muts.maf"),
                        "class": "File"
                    },
                    "snp_pileup": {
                         "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_DIR'], "Sample3.rg.md.abra.printreads__Sample4.rg.md.abra.printreads.dat.gz"),
                         "class": "File"
                    },
                    "pair_id": "Sample4.Sample3",
                    "tumor_id": "Sample4",
                    "normal_id": "Sample3"
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
                'analysis_dir': {
                    'class': 'Directory',
                    'basename': 'analysis',
                    'location': 'file://' + os.path.join(output_dir, 'analysis'),
                    'path': os.path.join(output_dir, 'analysis'),
                    'listing': [
                        {'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.gene.cna.txt'),
                        'basename': 'Proj_08390_G.gene.cna.txt',
                        'class': 'File',
                        'checksum': 'sha1$adf367afca81f3bb6e75272c86c38c4ad20246a9',
                        'size': 170207,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.gene.cna.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf'),
                        'basename': 'Proj_08390_G.muts.maf',
                        'class': 'File',
                        'checksum': 'sha1$38ccb013d912dbf8ea4d508ac112513a6e44c88c',
                        'size': 61631,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf')},
                        {'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.seg.cna.txt'),
                        'basename': 'Proj_08390_G.seg.cna.txt',
                        'class': 'File',
                        'checksum': 'sha1$8995133edb5b4547371b0fdea9925848eb480fbc',
                        'size': 2703,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.seg.cna.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.svs.maf'),
                        'basename': 'Proj_08390_G.svs.maf',
                        'class': 'File',
                        'checksum': 'sha1$5c2a63fc01980550108e58079a8b689d53c97d8c',
                        'size': 35595,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.svs.maf')}
                    ]
                },
                'facets_dir': {'basename': 'facets',
                  'class': 'Directory',
                  'listing': [{'basename': 'Sample1.Sample2',
                               'class': 'Directory',
                               'listing': [{'basename': 'Sample1.Sample2_hisens.ccf.portal.maf',
                                            'checksum': 'sha1$fa311e483303fbda7778622af7be57c7d62aac76',
                                            'class': 'File',
                                            'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.Sample2/Sample1.Sample2_hisens.ccf.portal.maf'),
                                            'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.Sample2/Sample1.Sample2_hisens.ccf.portal.maf'),
                                            'size': 19822273},
                                           {'basename': 'Sample1.arm_level.txt',
                                            'checksum': 'sha1$4f1d6620f3c613e40a147338739d3267e891dda7',
                                            'class': 'File',
                                            'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.arm_level.txt'),
                                            'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.arm_level.txt'),
                                            'size': 1986},
                                           {'basename': 'Sample1.txt',
                                            'checksum': 'sha1$8abe9bed5d1cb807b451c0647f658394d3f96708',
                                            'class': 'File',
                                            'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.txt'),
                                            'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.txt'),
                                            'size': 560},
                                           {'basename': 'Sample1.gene_level.txt',
                                            'checksum': 'sha1$6b31bff7758586ec39c94a0c2ddb8201eda5e0ce',
                                            'class': 'File',
                                            'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.gene_level.txt'),
                                            'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.gene_level.txt'),
                                            'size': 252793},
                                           {'basename': 'Sample1_hisens.cncf.txt',
                                            'checksum': 'sha1$a183b43abb9b7bb8f67ec6a0411539643d1275fe',
                                            'class': 'File',
                                            'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.cncf.txt'),
                                            'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.cncf.txt'),
                                            'size': 5774},
                                           {'basename': 'Sample1_hisens.rds',
                                            'checksum': 'sha1$d46031ce35efd7341146cac0db809c8cb75c97c2',
                                            'class': 'File',
                                            'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.rds'),
                                            'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.rds'),
                                            'size': 550930},
                                           {'basename': 'Sample1_hisens.seg',
                                            'checksum': 'sha1$1e0f1ad1cc751d53f02ad19717921027c3c96dfb',
                                            'class': 'File',
                                            'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.seg'),
                                            'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_hisens.seg'),
                                            'size': 1786},
                                           {'basename': 'Sample1_purity.rds',
                                            'checksum': 'sha1$265c7b94df391ae12e3bc1cc3ec4e763629a4f0a',
                                            'class': 'File',
                                            'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_purity.rds'),
                                            'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_purity.rds'),
                                            'size': 550388},
                                           {'basename': 'Sample1_purity.seg',
                                            'checksum': 'sha1$36904773ac5f194eac7f816a5fa9b055931bb48c',
                                            'class': 'File',
                                            'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_purity.seg'),
                                            'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1_purity.seg'),
                                            'size': 1284},
                                           {'basename': 'Sample1.qc.txt',
                                            'checksum': 'sha1$6fdefe96a23c0bc52050c46a16625fca603e4eb6',
                                            'class': 'File',
                                            'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.qc.txt'),
                                            'path': os.path.join(output_dir,'facets/Sample1.Sample2/Sample1.qc.txt'),
                                            'size': 1340}],
                               'location': 'file://' + os.path.join(output_dir,'facets/Sample1.Sample2'),
                               'path': os.path.join(output_dir,'facets/Sample1.Sample2')},
                              {'basename': 'Sample4.Sample3',
                               'class': 'Directory',
                               'listing': [{'basename': 'Sample4.Sample3_hisens.ccf.portal.maf',
                                            'checksum': 'sha1$86af92550b3e4252f950223f997579bfe651b65f',
                                            'class': 'File',
                                            'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3/Sample4.Sample3_hisens.ccf.portal.maf'),
                                            'path': os.path.join(output_dir,'facets/Sample4.Sample3/Sample4.Sample3_hisens.ccf.portal.maf'),
                                            'size': 1136054},
                                           {'basename': 'Sample4.arm_level.txt',
                                            'checksum': 'sha1$09b97c42a9875d642dae4d030712264f1b92b550',
                                            'class': 'File',
                                            'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3/Sample4.arm_level.txt'),
                                            'path': os.path.join(output_dir,'facets/Sample4.Sample3/Sample4.arm_level.txt'),
                                            'size': 1424},
                                           {'basename': 'Sample4.txt',
                                            'checksum': 'sha1$1bb6cf1ea21bc774037913633e5e06772d791b35',
                                            'class': 'File',
                                            'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3/Sample4.txt'),
                                            'path': os.path.join(output_dir,'facets/Sample4.Sample3/Sample4.txt'),
                                            'size': 602},
                                           {'basename': 'Sample4.gene_level.txt',
                                            'checksum': 'sha1$4a1f580aebf86a3a2ce4d5d27bc1caa3655e8788',
                                            'class': 'File',
                                            'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3/Sample4.gene_level.txt'),
                                            'path': os.path.join(output_dir,'facets/Sample4.Sample3/Sample4.gene_level.txt'),
                                            'size': 144319},
                                           {'basename': 'Sample4_hisens.cncf.txt',
                                            'checksum': 'sha1$a51eac1442aac591412b79305f003d464322140e',
                                            'class': 'File',
                                            'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3/Sample4_hisens.cncf.txt'),
                                            'path': os.path.join(output_dir,'facets/Sample4.Sample3/Sample4_hisens.cncf.txt'),
                                            'size': 4802},
                                           {'basename': 'Sample4_hisens.rds',
                                            'checksum': 'sha1$33286048372616b019fcf52587b25bbbf5d1ae08',
                                            'class': 'File',
                                            'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3/Sample4_hisens.rds'),
                                            'path': os.path.join(output_dir,'facets/Sample4.Sample3/Sample4_hisens.rds'),
                                            'size': 186732},
                                           {'basename': 'Sample4_hisens.seg',
                                            'checksum': 'sha1$ea32de058625243fdc9f48d9cbab746444549e2e',
                                            'class': 'File',
                                            'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3/Sample4_hisens.seg'),
                                            'path': os.path.join(output_dir,'facets/Sample4.Sample3/Sample4_hisens.seg'),
                                            'size': 1584},
                                           {'basename': 'Sample4_purity.rds',
                                            'checksum': 'sha1$1eda8dc016cf8322cab1190ea4f723ecb33c4440',
                                            'class': 'File',
                                            'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3/Sample4_purity.rds'),
                                            'path': os.path.join(output_dir,'facets/Sample4.Sample3/Sample4_purity.rds'),
                                            'size': 186513},
                                           {'basename': 'Sample4_purity.seg',
                                            'checksum': 'sha1$6be6775133173d609e176f16a70c5a3c2aeb9027',
                                            'class': 'File',
                                            'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3/Sample4_purity.seg'),
                                            'path': os.path.join(output_dir,'facets/Sample4.Sample3/Sample4_purity.seg'),
                                            'size': 1276},
                                           {'basename': 'Sample4.qc.txt',
                                            'checksum': 'sha1$0c17bd63f43a6726682dec93a58cbcb3ccfd20e1',
                                            'class': 'File',
                                            'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3/Sample4.qc.txt'),
                                            'path': os.path.join(output_dir,'facets/Sample4.Sample3/Sample4.qc.txt'),
                                            'size': 1330}],
                               'location': 'file://' + os.path.join(output_dir,'facets/Sample4.Sample3'),
                               'path': os.path.join(output_dir,'facets/Sample4.Sample3')}],
                  'location': 'file://' + os.path.join(output_dir,'facets'),
                  'path': os.path.join(output_dir,'facets')
                },
                'facets_failed_pairs': [],
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
                        'checksum': 'sha1$dd6652702fbe4b4390bc5d209ce6f4a4adf7e17a',
                        'size': 8560,
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
                        'checksum': 'sha1$711e5880562e5571d5c9d6538de59ef7b6e746c3',
                        'size': 6741,
                        'path': os.path.join(output_dir, 'portal/data_CNA.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/data_CNA.ascna.txt'),
                        'basename': 'data_CNA.ascna.txt',
                        'class': 'File',
                        'checksum': 'sha1$74269b7694b9c5572a5a7e94411756e0418665c0',
                        'size': 8810,
                        'path': os.path.join(output_dir, 'portal/data_CNA.ascna.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/data_mutations_extended.txt'),
                        'basename': 'data_mutations_extended.txt',
                        'class': 'File',
                        'checksum': 'sha1$8e12ceda93a0b43c1a98c0903269b73dff7438c8',
                        'size': 8344,
                        'path': os.path.join(output_dir, 'portal/data_mutations_extended.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'portal/Proj_08390_G_data_cna_hg19.seg'),
                        'basename': 'Proj_08390_G_data_cna_hg19.seg',
                        'class': 'File',
                        'checksum': 'sha1$8995133edb5b4547371b0fdea9925848eb480fbc',
                        'size': 2703,
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
            comments, mutations = load_mutations(os.path.join(output_dir, 'analysis', 'Proj_08390_G.muts.maf'))
            self.assertEqual(len(mutations), 34)
            comments, mutations = load_mutations(os.path.join(output_dir, 'portal', 'data_mutations_extended.txt'))
            self.assertEqual(len(mutations), 27)

if __name__ == "__main__":
    unittest.main()
