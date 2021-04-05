#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for running a CWL with a pre-made JSON input file
"""
import argparse
import glob
import os
from pluto.settings import CWL_DIR
from pluto.tools import CWLRunner

# get list of all CWL files in the cwl dir
cwl_filepaths = glob.glob(CWL_DIR + '/*.cwl')
cwl_files = { os.path.basename(filepath):filepath for filepath in cwl_filepaths }

def main(**kwargs):
    """
    Run the CWL from a JSON input
    """
    kwargs.pop('func') # comes from argparse
    cwl_file = kwargs.pop('cwl_file')
    cwl_filepath = cwl_files[cwl_file]
    input_json = kwargs.pop('input_json')
    # print(cwl_file, input_json, kwargs)
    runner = CWLRunner(cwl_file = cwl_filepath, input = input_json, input_is_file = True, **kwargs)
    output_json, output_dir, output_json_file = runner.run()
    return(output_json, output_dir, output_json_file)

def parse():
    """
    CLI argument parsing
    """
    parser = argparse.ArgumentParser(description = '')
    parser.add_argument("--engine", default = 'cwltool', dest = 'engine', choices = ['cwltool', 'toil'], help = "CWL execution engine to use")
    parser.add_argument("--print-command", action = 'store_true', dest = 'print_command', help = "Print the CWL runner command and exit")
    parser.add_argument("--restart", action = 'store_true', dest = 'restart', help = "Restart a previous run; requires jobStore")
    parser.add_argument("--debug", action = 'store_true', dest = 'debug', help = "Restart a previous run; requires jobStore")
    parser.add_argument("--jobStore", dest = 'jobStore', default = None, help = "Job store to use for a restarted run")
    parser.add_argument("--dir", dest = 'dir', default = None, help = "Directory where the CWL will be executed")
    parser.add_argument("--output-dir", dest = 'output_dir', default = None, help = "Directory where the CWL output will be saved")
    parser.add_argument("--parallel", dest = 'parallel', action = "store_true", help = "Run cwltool in parallel mode. NOTE: make sure all containers are cached")
    parser.add_argument("--cwl", dest = 'cwl_file', required = True, choices = cwl_files.keys(), help = "CWL filename to be executed")
    parser.add_argument('input_json', help="Input JSON file to run the CWL with")

    parser.set_defaults(func = main)

    args = parser.parse_args()
    args.func(**vars(args))

if __name__ == '__main__':
    parse()
