import json
import copy
from pluto.tools import CWLRunner

class Operator(object):
    """
    Object class to handle argument parsing for running CWL's
    Meant to be subclassed, then imported into run.py and invoked from the self._run() entrypoint
    All required CWL inputs that do not need special processing should be passed in from the run.main argument parser;
    they will be loaded into Operator.args, then they should be converted into the correct CWL format inside self.__init__
    """
    # placeholder values to be overriden by subclasses
    cwl_file = None # str or CWLFile object
    input = {}
    pair_template = {}
    print_input = False
    runner_args = dict(
        engine = 'cwltool',
        verbose = True,
        dir = None,
        output_dir = None,
        print_command = False,
        restart = False,
        jobStore = None,
        debug = False
        )

    def __init__(self, **kwargs):
        self.args = copy.deepcopy(kwargs)
        # use these keyword args to set instance attributes
        if 'dir' in self.args:
            self.runner_args['dir'] = self.args.pop('dir')
        if 'output_dir' in self.args:
            self.runner_args['output_dir'] = self.args.pop('output_dir')
        if 'verbose' in self.args:
            self.runner_args['verbose'] = self.args.pop('verbose')
        if 'print_command' in self.args:
            self.runner_args['print_command'] = self.args.pop('print_command')
        if 'engine' in self.args:
            self.runner_args['engine'] = self.args.pop('engine')
        if 'restart' in self.args:
            self.runner_args['restart'] = self.args.pop('restart')
        if 'jobStore' in self.args:
            self.runner_args['jobStore'] = self.args.pop('jobStore')
        if 'debug' in self.args:
            self.runner_args['debug'] = self.args.pop('debug')



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

        runner = CWLRunner(cwl_file = self.cwl_file, input = self.input, **self.runner_args)
        output_json, output_dir, output_json_file = runner.run()
        return(output_json, output_dir, output_json_file)

    def __str__(self):
        return(str(vars(self)))
