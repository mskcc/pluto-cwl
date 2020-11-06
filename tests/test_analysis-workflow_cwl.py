#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unit tests for the workflow.cwl
"""
import os
import json
import unittest
from tempfile import TemporaryDirectory, NamedTemporaryFile


# relative imports, from CLI and from parent project
if __name__ != "__main__":
    from .tools import run_cwl, load_mutations
    from .settings import CWL_DIR, DATA_SETS

if __name__ == "__main__":
    from tools import run_cwl, load_mutations
    from settings import CWL_DIR, DATA_SETS

cwl_file = os.path.join(CWL_DIR, 'analysis-workflow.cwl')

class TestAnalysisWorkflow(unittest.TestCase):
    def test_run_worflow_one_maf(self):
        """
        Test that the workflow works correctly when run with a single maf
        """
        input_json = {
            "is_impact": True,
            "argos_version_string": "2.x",
            "analysis_gene_cna_filename": "Proj_08390_G.gene.cna.txt",
            "analysis_mutations_filename": "Proj_08390_G.muts.maf",
            "analysis_mutations_share_filename": "Proj_08390_G.muts.share.maf",
            "analysis_segment_cna_filename": "Proj_08390_G.seg.cna.txt",
            "analysis_sv_filename": "Proj_08390_G.svs.maf",
            "helix_filter_version": "20.06.1",
            "targets_list": {
                "path": DATA_SETS['Proj_08390_G']["targets_list"],
                "class": "File"
            },
            "mutation_maf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf"),
                    "class": "File"
                }
            ],
            "mutation_svs_maf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.svs.pass.vep.maf"),
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
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file)

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
                            'checksum': 'sha1$7cc89d24556de93b9a409812317581e67e5df494',
                            'size': 87905,
                            'path': os.path.join(output_dir, 'analysis/Proj_08390_G.gene.cna.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf'),
                            'basename': 'Proj_08390_G.muts.maf',
                            'class': 'File',
                            'checksum': 'sha1$25247bee2d3fdf76dd64426e327736798660c6b4',
                            'size': 32346,
                            'path': os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.muts.share.maf'),
                            'basename': 'Proj_08390_G.muts.share.maf',
                            'class': 'File',
                            'checksum': 'sha1$86a4932b612b73752ecdcda111d32b91f6c937d9',
                            'size': 6565,
                            'path': os.path.join(output_dir, 'analysis/Proj_08390_G.muts.share.maf')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.seg.cna.txt'),
                            'basename': 'Proj_08390_G.seg.cna.txt',
                            'class': 'File',
                            'checksum': 'sha1$f0ebb82c34b6530447fa1e70b6dedcc039840d61',
                            'size': 1632,
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
                    }
                }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)
            comments, mutations = load_mutations(os.path.join(output_dir, 'analysis', 'Proj_08390_G.muts.maf'))
            self.assertEqual(len(mutations), 22)

    def test_run_worflow_two_mafs(self):
        """
        Test that the workflow works correctly when run with two maf files
        """
        input_json = {
            "is_impact": True,
            "argos_version_string": "2.x",
            "analysis_gene_cna_filename": "Proj_08390_G.gene.cna.txt",
            "analysis_mutations_filename": "Proj_08390_G.muts.maf",
            "analysis_mutations_share_filename": "Proj_08390_G.muts.share.maf",
            "analysis_segment_cna_filename": "Proj_08390_G.seg.cna.txt",
            "analysis_sv_filename": "Proj_08390_G.svs.maf",
            "helix_filter_version": "20.06.1",
            "targets_list": {
                "path": DATA_SETS['Proj_08390_G']["targets_list"],
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
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file)

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
                        'checksum': 'sha1$ab17d587ad5ae0a87fd6c6d4dd2d5d1701208ce9',
                        'size': 173982,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.gene.cna.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf'),
                        'basename': 'Proj_08390_G.muts.maf',
                        'class': 'File',
                        'checksum': 'sha1$61a99e4c9ca4dbce12e52a3ed1c635738a162ffd',
                        'size': 53071,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf')},
                        {
                            'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.muts.share.maf'),
                            'basename': 'Proj_08390_G.muts.share.maf',
                            'class': 'File',
                            'checksum': 'sha1$c579fa93d8b096567d5d34bab09db9a5e1757302',
                            'size': 9569,
                            'path': os.path.join(output_dir, 'analysis/Proj_08390_G.muts.share.maf')
                        },
                        {'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.seg.cna.txt'),
                        'basename': 'Proj_08390_G.seg.cna.txt',
                        'class': 'File',
                        'checksum': 'sha1$f6a77b280c047a7e2082e3a09e8138f861790d3a',
                        'size': 3191,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.seg.cna.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.svs.maf'),
                        'basename': 'Proj_08390_G.svs.maf',
                        'class': 'File',
                        'checksum': 'sha1$5c2a63fc01980550108e58079a8b689d53c97d8c',
                        'size': 35595,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.svs.maf')}
                    ]
                }
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)
            comments, mutations = load_mutations(os.path.join(output_dir, 'analysis', 'Proj_08390_G.muts.maf'))
            self.assertEqual(len(mutations), 34)

    def test_run_worflow_mixed_mafs(self):
        """
        Test that the workflow works correctly when run with a mix of Argos muts.maf files and Facets Suite annotated maf files
        The Facets Suite maf files have extra columns that need to be retained in the output
        """
        input_json = {
            "is_impact": True,
            "argos_version_string": "2.x",
            "analysis_gene_cna_filename": "Proj_08390_G.gene.cna.txt",
            "analysis_mutations_filename": "Proj_08390_G.muts.maf",
            "analysis_mutations_share_filename": "Proj_08390_G.muts.share.maf",
            "analysis_segment_cna_filename": "Proj_08390_G.seg.cna.txt",
            "analysis_sv_filename": "Proj_08390_G.svs.maf",
            "helix_filter_version": "20.06.1",
            "targets_list": {
                "path": DATA_SETS['Proj_08390_G']["targets_list"],
                "class": "File"
            },
            "mutation_maf_files": [
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['MAF_DIR'], "Sample1.Sample2.muts.maf"),
                    "class": "File"
                },
                {
                    "path": os.path.join(DATA_SETS['Proj_08390_G']['FACETS_SUITE_DIR'], "Sample4.Sample3_hisens.ccf.portal.maf"),
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
            output_json, output_dir = run_cwl(
                testcase = self,
                tmpdir = tmpdir,
                input_json = input_json,
                cwl_file = cwl_file)


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
                        'checksum': 'sha1$ab17d587ad5ae0a87fd6c6d4dd2d5d1701208ce9',
                        'size': 173982,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.gene.cna.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf'),
                        'basename': 'Proj_08390_G.muts.maf',
                        'class': 'File',
                        'checksum': 'sha1$462310d9ca98e9475b75e601c75d2c6e2c3cf2ec',
                        'size': 58528,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf')},
                        {
                            'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.muts.share.maf'),
                            'basename': 'Proj_08390_G.muts.share.maf',
                            'class': 'File',
                            'checksum': 'sha1$5770ed1b1e386ce30b210b30c6cd3918868cf299',
                            'size': 9342,
                            'path': os.path.join(output_dir, 'analysis/Proj_08390_G.muts.share.maf')
                        },
                        {'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.seg.cna.txt'),
                        'basename': 'Proj_08390_G.seg.cna.txt',
                        'class': 'File',
                        'checksum': 'sha1$f6a77b280c047a7e2082e3a09e8138f861790d3a',
                        'size': 3191,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.seg.cna.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.svs.maf'),
                        'basename': 'Proj_08390_G.svs.maf',
                        'class': 'File',
                        'checksum': 'sha1$5c2a63fc01980550108e58079a8b689d53c97d8c',
                        'size': 35595,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.svs.maf')}
                    ]
                }
            }
            self.maxDiff = None
            self.assertDictEqual(output_json, expected_output)

            comments, mutations = load_mutations(os.path.join(output_dir, 'analysis', 'Proj_08390_G.muts.maf'))
            self.assertEqual(len(mutations), 34)

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

            self.assertEqual(mutations[0]['t_af'], '0.42953020134228187')

if __name__ == "__main__":
    unittest.main()
