#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Operator to run the workflow with Facets based on CLI input
"""
import os
import sys
import json
import subprocess
from pluto.tools import CWLFile, CWLRunner
from pluto.settings import DATA_SETS, KNOWN_FUSIONS_FILE, IMPACT_FILE
from .input import generate_input

cwl_file = CWLFile('workflow_with_facets.cwl')

# try to get the git branch and version
version = None
try:
    version = subprocess.check_output(['git', 'describe', '--tags']).decode('ascii').strip()
except:
    pass

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

def main(
    assay_coverage,
    project_id,
    cancer_type,
    pairs_file,
    data_clinical_file,
    sample_summary_file,
    print_input = False,
    dir = None,
    verbose = True,
    **kwargs):
    # collect and parse through all the input args
    args = {**kwargs}
    if 'func' in args:
        args.pop('func')
    # required args
    args['assay_coverage'] = assay_coverage
    args['project_id'] = project_id
    args['cancer_type'] = cancer_type
    # samplesheets
    args['pairs_file'] = os.path.abspath(pairs_file)
    args['data_clinical_file'] = os.path.abspath(data_clinical_file)
    args['sample_summary_file'] = os.path.abspath(sample_summary_file)
    # file list sheets
    args['mutation_svs_txt_files'] = args.get('mutation_svs_txt_files', None)
    args['mutation_svs_maf_files'] = args.get('mutation_svs_maf_files', None)

    # optional args; { 'arg': 'default'}
    opt_args = {
        'cancer_study_identifier': project_id,
        'project_name': project_id,
        'project_short_name': project_id,
        'project_description': '',
        'project_pi': '',
        'request_pi': '',
        'is_impact': True,
        'argos_version_string': '2.x',
        'analysis_gene_cna_filename': project_id + '.gene.cna.txt',
        'analysis_mutations_filename': project_id + '.muts.maf',
        'analysis_mutations_share_filename': project_id + '.muts.share.maf',
        'analysis_segment_cna_filename': project_id + '.seg.cna.txt',
        'analysis_sv_filename': project_id + '.svs.maf',
        'cbio_meta_cna_segments_filename': project_id + '_meta_cna_hg19_seg.txt',
        'cbio_segment_data_filename': project_id + '_data_cna_hg19.seg',
        'helix_filter_version': version,
        'IMPACT_gene_list': IMPACT_FILE,
        'targets_list': DATA_SETS['Proj_08390_G']["targets_list"],
        'known_fusions_file': KNOWN_FUSIONS_FILE
        }
    # update missing args with default value
    for key, value in opt_args.items():
        if key not in args:
            args[key] = value
        elif args[key] is None:
            args[key] = value

    input = generate_input(
        args,
        pair_template = pair_template,
        bool_keys = ['is_impact'],
        File_keys = ['IMPACT_gene_list', 'data_clinical_file', 'sample_summary_file', 'targets_list', 'known_fusions_file'],
        array_File_keys = ['mutation_svs_txt_files', 'mutation_svs_maf_files']
    )

    if print_input:
        print(json.dumps(input, indent = 4))
        return()

    runner = CWLRunner(cwl_file, input, dir = dir, verbose = verbose)
    runner.run()
