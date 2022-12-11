#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
unit tests for the workflow.cwl
"""
import os
import sys
import unittest

PARENT_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.insert(0, PARENT_DIR)
from pluto import (
    CWLFile, 
    PlutoTestCase,
    DATA_SETS,
    IMPACT_FILE,
    OFile,
    ODir
)
sys.path.pop(0)

class TestAnalysisWorkflow(PlutoTestCase):
    cwl_file = CWLFile('analysis-workflow.cwl')
    def test_run_worflow_one_maf(self):
        """
        Test that the workflow works correctly when run with a single maf
        """
        self.input = {
            "is_impact": True,
            "argos_version_string": "2.x",
            "analysis_gene_cna_filename": "Proj_08390_G.gene.cna.txt",
            "analysis_mutations_filename": "Proj_08390_G.muts.maf",
            "analysis_mutations_share_filename": "Proj_08390_G.muts.share.maf",
            "analysis_segment_cna_filename": "Proj_08390_G.seg.cna.txt",
            "analysis_sv_filename": "Proj_08390_G.svs.maf",
            "helix_filter_version": "20.06.1",
            "IMPACT_gene_list": {
                  "class": "File",
                  "path": IMPACT_FILE
                },
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

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'analysis_dir': ODir(name='analysis', items=[
                OFile(name='Proj_08390_G.gene.cna.txt', size=87905, hash='7cc89d24556de93b9a409812317581e67e5df494'),
                OFile(name='Proj_08390_G.muts.maf', size=33243, hash='2c8904927a917d6e935ef207582d995680574d16'),
                OFile(name='Proj_08390_G.muts.share.maf', size=7462, hash='b5af4e0fcd89fecabf8095aa3d7690e5edb8dca1'),
                OFile(name='Proj_08390_G.seg.cna.txt', size=1632, hash='f0ebb82c34b6530447fa1e70b6dedcc039840d61'),
                OFile(name='Proj_08390_G.svs.maf', size=23603, hash='df420706bb5b772a79317843c0a01a3c88a9571d')], dir=output_dir)
            }
        self.assertCWLDictEqual(output_json, expected_output)
        self.assertNumMutations(os.path.join(output_dir, 'analysis', 'Proj_08390_G.muts.maf'), 22)

    def test_run_worflow_two_mafs(self):
        """
        Test that the workflow works correctly when run with two maf files
        """
        self.input = {
            "is_impact": True,
            "argos_version_string": "2.x",
            "analysis_gene_cna_filename": "Proj_08390_G.gene.cna.txt",
            "analysis_mutations_filename": "Proj_08390_G.muts.maf",
            "analysis_mutations_share_filename": "Proj_08390_G.muts.share.maf",
            "analysis_segment_cna_filename": "Proj_08390_G.seg.cna.txt",
            "analysis_sv_filename": "Proj_08390_G.svs.maf",
            "helix_filter_version": "20.06.1",
            "IMPACT_gene_list": {
                  "class": "File",
                  "path": IMPACT_FILE
                },
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

        output_json, output_dir = self.run_cwl()

        expected_output = {
            'analysis_dir': ODir(name='analysis', items=[
                OFile(name='Proj_08390_G.gene.cna.txt', size=173982, hash='ab17d587ad5ae0a87fd6c6d4dd2d5d1701208ce9'),
                OFile(name='Proj_08390_G.muts.maf', size=54458, hash='d4352ee2b702877b84db2b632972ccad2441f3e0'),
                OFile(name='Proj_08390_G.muts.share.maf', size=10956, hash='086ce6517eae68e47160c8740c5f00d7c3454110'),
                OFile(name='Proj_08390_G.seg.cna.txt', size=3191, hash='f6a77b280c047a7e2082e3a09e8138f861790d3a'),
                OFile(name='Proj_08390_G.svs.maf', size=35595, hash='5c2a63fc01980550108e58079a8b689d53c97d8c')], dir=output_dir)
        }

        self.assertCWLDictEqual(output_json, expected_output)
        self.assertNumMutations(os.path.join(output_dir, 'analysis', 'Proj_08390_G.muts.maf'), 34)

    def test_run_worflow_mixed_mafs(self):
        """
        Test that the workflow works correctly when run with a mix of Argos muts.maf files and Facets Suite annotated maf files
        The Facets Suite maf files have extra columns that need to be retained in the output
        """
        self.input = {
            "is_impact": True,
            "argos_version_string": "2.x",
            "analysis_gene_cna_filename": "Proj_08390_G.gene.cna.txt",
            "analysis_mutations_filename": "Proj_08390_G.muts.maf",
            "analysis_mutations_share_filename": "Proj_08390_G.muts.share.maf",
            "analysis_segment_cna_filename": "Proj_08390_G.seg.cna.txt",
            "analysis_sv_filename": "Proj_08390_G.svs.maf",
            "helix_filter_version": "20.06.1",
            "IMPACT_gene_list": {
                  "class": "File",
                  "path": IMPACT_FILE
                },
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

        output_json, output_dir = self.run_cwl()

        expected_output = {
        'analysis_dir': ODir(name='analysis', items=[
            OFile(name='Proj_08390_G.gene.cna.txt', size=173982, hash='ab17d587ad5ae0a87fd6c6d4dd2d5d1701208ce9'),
            OFile(name='Proj_08390_G.muts.maf', size=59915, hash='66a87cb8cc2eea31f490852d468bedd958c4ecc5'),
            OFile(name='Proj_08390_G.muts.share.maf', size=10729, hash='cbaa23bb848978cde135efd3870db8f35b3f2861'),
            OFile(name='Proj_08390_G.seg.cna.txt', size=3191, hash='f6a77b280c047a7e2082e3a09e8138f861790d3a'),
            OFile(name='Proj_08390_G.svs.maf', size=35595, hash='5c2a63fc01980550108e58079a8b689d53c97d8c')], dir=output_dir)
        }

        self.assertCWLDictEqual(output_json, expected_output)

        self.assertNumMutations(os.path.join(output_dir, 'analysis', 'Proj_08390_G.muts.maf'), 34)


        # colnames = mutations[0].keys()
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
        self.assertMutHeadersContain(os.path.join(output_dir, 'analysis', 'Proj_08390_G.muts.maf'), some_required_colnames)

        # TODO: reimplement these tests;
        # self.assertEqual(mutations[0]['t_af'], '0.42953020134228187')
        # self.assertEqual(mutations[0]['is_in_impact'], 'True')
        # self.assertEqual(mutations[0]['impact_assays'], 'IMPACT341,IMPACT410,IMPACT468,IMPACT505')




