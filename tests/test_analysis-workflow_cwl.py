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
    from .tools import run_command, load_mutations
    from .settings import CWL_DIR, CWL_ARGS, DATA_SETS, KNOWN_FUSIONS_FILE

if __name__ == "__main__":
    from tools import run_command, load_mutations
    from settings import CWL_DIR, CWL_ARGS, DATA_SETS, KNOWN_FUSIONS_FILE

cwl_file = os.path.join(CWL_DIR, 'analysis-workflow.cwl')

class TestAnalysisWorkflow(unittest.TestCase):
    def test_run_worflow_one_maf(self):
        """
        Test that the workflow works correctly when run with a single maf
        """
        input_json = {
            "is_impact": "True",
            "argos_version_string": "2.x",
            "analysis_gene_cna_filename": "Proj_08390_G.gene.cna.txt",
            "analysis_mutations_filename": "Proj_08390_G.muts.maf",
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
                            'checksum': 'sha1$7cc89d24556de93b9a409812317581e67e5df494',
                            'size': 87905,
                            'path': os.path.join(output_dir, 'analysis/Proj_08390_G.gene.cna.txt')
                        },
                        {
                            'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf'),
                            'basename': 'Proj_08390_G.muts.maf',
                            'class': 'File',
                            'checksum': 'sha1$9b506f93b744ea6dd35c7ef9ba29b4f91cb92de1',
                            'size': 31943,
                            'path': os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf')
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

    def test_run_worflow_two_mafs(self):
        """
        Test that the workflow works correctly when run with two maf files
        """
        input_json = {
            "is_impact": "True",
            "argos_version_string": "2.x",
            "analysis_gene_cna_filename": "Proj_08390_G.gene.cna.txt",
            "analysis_mutations_filename": "Proj_08390_G.muts.maf",
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
                        'checksum': 'sha1$ab17d587ad5ae0a87fd6c6d4dd2d5d1701208ce9',
                        'size': 173982,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.gene.cna.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf'),
                        'basename': 'Proj_08390_G.muts.maf',
                        'class': 'File',
                        'checksum': 'sha1$aea0871b7acf07ced0621340c78a41bad52c20bc',
                        'size': 52475,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf')},
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

    def test_run_worflow_mixed_mafs(self):
        """
        Test that the workflow works correctly when run with a mix of Argos muts.maf files and Facets Suite annotated maf files
        The Facets Suite maf files have extra columns that need to be retained in the output
        """
        input_json = {
            "is_impact": "True",
            "argos_version_string": "2.x",
            "analysis_gene_cna_filename": "Proj_08390_G.gene.cna.txt",
            "analysis_mutations_filename": "Proj_08390_G.muts.maf",
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
                        'checksum': 'sha1$ab17d587ad5ae0a87fd6c6d4dd2d5d1701208ce9',
                        'size': 173982,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.gene.cna.txt')},
                        {'location': 'file://' + os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf'),
                        'basename': 'Proj_08390_G.muts.maf',
                        'class': 'File',
                        'checksum': 'sha1$1a964935bb062e9c70a8dd02ca49ad82a9278af2',
                        'size': 57919,
                        'path': os.path.join(output_dir, 'analysis/Proj_08390_G.muts.maf')},
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


if __name__ == "__main__":
    unittest.main()
