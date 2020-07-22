"""
Helper functions for running tests
"""
import subprocess as sp
import csv

def run_command(args):
    """
    Helper function to run a shell command easier

    Parameters
    ----------
    args: list
        a list of shell args to execute
    """
    process = sp.Popen(args, stdout = sp.PIPE, stderr = sp.PIPE, universal_newlines = True)
    proc_stdout, proc_stderr = process.communicate()
    returncode = process.returncode
    proc_stdout = proc_stdout.strip()
    proc_stderr = proc_stderr.strip()
    return(returncode, proc_stdout, proc_stderr)

def parse_header_comments(filename):
    """
    Parse a file with comments in its header to return the comments and the line number to start reader from

    comments, start_line = parse_header_comments(filename)
    with open(portal_file) as fin:
        while start_line > 0:
            next(fin)
            start_line -= 1
        reader = csv.DictReader(fin, delimiter = '\t') # header_line = next(fin)
        portal_lines = [ row for row in reader ]
    """
    comments = []
    start_line = 0
    # find the first line without comments
    with open(filename) as fin:
        for i, line in enumerate(fin):
            if line.startswith('#'):
                comments.append(line.strip())
                start_line += 1
    return(comments, start_line)

def load_mutations(filename):
    """
    Load the mutations from a file to use for testing
    """
    comments, start_line = parse_header_comments(filename)
    with open(filename) as fin:
        while start_line > 0:
            next(fin)
            start_line -= 1
        reader = csv.DictReader(fin, delimiter = '\t')
        mutations = [ row for row in reader ]
    return(comments, mutations)
