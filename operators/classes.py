from pluto.tools import CWLRunner

class Operator(object):
    # placeholder values to be overriden by subclasses
    cwl_file = None # str or CWLFile object
    dir = 'pipeline_output'
    verbose = True
    input = None

    def __init__(self, **kwargs):
        pass
        # for k,v in kwargs.items():
        #     if k not in ['func']:
        #         setattr(self, k, v)

    def parse_kwargs(self, **kwargs):
        """
        Override this method to parse the args passed in a set instance attributes
        """
        if 'dir' in kwargs:
            self.dir = kwargs['dir']
        if 'verbose' in kwargs:
            self.verbose = kwargs['verbose']

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
        runner = CWLRunner(
            cwl_file = self.cwl_file,
            input = self.input,
            dir = self.dir,
            verbose = self.verbose)
        output_json, output_dir, output_json_file = runner.run()
        return(output_json, output_dir, output_json_file)

    def __str__(self):
        return(str(vars(self)))
