#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run the workflow with Facets based on CLI input
"""
import os
from pluto.tools import CWLFile
from pluto.settings import DATA_SETS, KNOWN_FUSIONS_FILE, IMPACT_FILE
from .input import generate_input
from .classes import Operator
import subprocess
import copy

# try to get the git branch and version
version = None
try:
    version = subprocess.check_output(['git', 'describe', '--tags']).decode('ascii').strip()
except:
    pass

class WorkflowWithFacets(Operator):
    cwl_file = CWLFile('workflow_with_facets.cwl')
    pair_template = {
        "pair_maf": {
            "path": None,
            "class": "File"
        },
        "snp_pileup": {
            "path": None,
            "class": "File"
        },
        "pair_id": None,
        "tumor_id": None,
        "normal_id": None
    }
    
    def __init__(self,
        assay_coverage,
        project_id,
        cancer_type,
        pairs_file,
        data_clinical_file,
        sample_summary_file,
        mutation_svs_txt_files = None,
        mutation_svs_maf_files = None,
        IMPACT_gene_list = IMPACT_FILE,
        targets_list = DATA_SETS['Proj_08390_G']["targets_list"],
        known_fusions_file = KNOWN_FUSIONS_FILE,
        version = version,
        **kwargs):
        super().__init__(**kwargs)
        # handling for required values with defaults
        self.args['pairs_file'] = os.path.abspath(pairs_file)
        self.args['data_clinical_file'] = os.path.abspath(data_clinical_file)
        self.args['sample_summary_file'] = os.path.abspath(sample_summary_file)
        self.args['project_id'] = project_id
        self.args['assay_coverage'] = assay_coverage
        self.args['cancer_type'] = cancer_type
        self.args['mutation_svs_txt_files'] = mutation_svs_txt_files
        self.args['mutation_svs_maf_files'] = mutation_svs_maf_files
        self.args['IMPACT_gene_list'] = IMPACT_gene_list
        self.args['targets_list'] = targets_list
        self.args['known_fusions_file'] = known_fusions_file
        self.version = version
        self.generate_input_data()

    def generate_input_data(self): # **kwargs
        # collect and parse through all the input args
        input_args = copy.deepcopy(self.args)

        # handling of optional args with defaults
        opt_args = { # 'arg': 'default'
            'cancer_study_identifier': self.args['project_id'],
            'project_name': self.args['project_id'],
            'project_short_name': self.args['project_id'],
            'project_description': '',
            'project_pi': '',
            'request_pi': '',
            'is_impact': True,
            'argos_version_string': '2.x',
            'analysis_gene_cna_filename': self.args['project_id'] + '.gene.cna.txt',
            'analysis_mutations_filename': self.args['project_id'] + '.muts.maf',
            'analysis_mutations_share_filename': self.args['project_id'] + '.muts.share.maf',
            'analysis_segment_cna_filename': self.args['project_id'] + '.seg.cna.txt',
            'analysis_sv_filename': self.args['project_id'] + '.svs.maf',
            'cbio_meta_cna_segments_filename': self.args['project_id'] + '_meta_cna_hg19_seg.txt',
            'cbio_segment_data_filename': self.args['project_id'] + '_data_cna_hg19.seg',
            'helix_filter_version': self.version
            }

        # update missing args with default value
        for key, value in opt_args.items():
            if key not in input_args:
                input_args[key] = value
            elif input_args[key] is None:
                input_args[key] = value

        self.input = generate_input(
            input_args,
            pair_template = self.pair_template,
            bool_keys = ['is_impact'],
            File_keys = ['IMPACT_gene_list', 'data_clinical_file', 'sample_summary_file', 'targets_list', 'known_fusions_file'],
            array_File_keys = ['mutation_svs_txt_files', 'mutation_svs_maf_files']
        )
