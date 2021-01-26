import json
import copy
from pluto.tools import CWLRunner

class Operator(object):
    # placeholder values to be overriden by subclasses
    cwl_file = None # str or CWLFile object
    dir = 'pipeline_output'
    verbose = True
    input = None
    print_input = False

    def __init__(self, **kwargs):
        self.args = copy.deepcopy(kwargs)
        # use these keyword args to set instance attributes
        if 'dir' in self.args:
            self.dir = self.args.pop('dir')
        if 'verbose' in self.args:
            self.verbose = self.args.pop('verbose')
        if 'print_input' in self.args:
            self.print_input = self.args.pop('print_input')
        if 'func' in self.args:
            self.args.pop('func')
        # pass
        # for k,v in kwargs.items():
        #     if k not in ['func']:
        #         setattr(self, k, v)

    def parse_kwargs(self, **kwargs):
        """
        Pre-processing of passed keyword args to set required instance values and strip out values that should not get passed to pipeline input later
        """

    def generate_input_data(self):
        """
        Override this method to turn instance attributes into the data for input to the pipeline
        """

    @classmethod
    def _run(cls, **kwargs):
        """
        Run entrypoint for CLI arg parser
        """
        instance = cls(**kwargs)
        instance.run()

    def run(self):
        """
        Run the CWL
        """
        if self.print_input:
            print(json.dumps(self.input, indent = 4))
            return()

        runner = CWLRunner(
            cwl_file = self.cwl_file,
            input = self.input,
            dir = self.dir,
            verbose = self.verbose)
        output_json, output_dir, output_json_file = runner.run()
        return(output_json, output_dir, output_json_file)

    def __str__(self):
        return(str(vars(self)))
