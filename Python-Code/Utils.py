

from timeit import default_timer
import argparse
import os


# --------------------------------------------------------------------------------------------------

class timer:

    # ----------------------------------------------------------------------------------------------

    def __init__(self,label,units='s'):
        """
        small tool for timing and printing timing info
        """

        self.label = label
        if units == 'm':
            self.units = 'm'
            self.scale = 1/60
        else:
            self.units = 's'
            self.scale = 1
        self.start_time = default_timer()

    # ----------------------------------------------------------------------------------------------

    def stop(self):
        """
        stop timer and print timing info
        """

        elapsed_time = default_timer()-self.start_time
        elapsed_time *= self.scale
        msg = '\n---------------------------------------------------------------------------\n'
        msg += f'timing:   {self.label} {elapsed_time:9.5f} [{self.units}]\n'
        print(msg)

# --------------------------------------------------------------------------------------------------




# --------------------------------------------------------------------------------------------------

class arg_parser:

    # ----------------------------------------------------------------------------------------------

    def __init__(self):
        """
        parses command line args
        """

        description = 'command line args for phonon explorer'
        cmd_parser = argparse.ArgumentParser(description=description)

        # input file
        help_msg = 'input file for phonon explorer'
        cmd_parser.add_argument('-i','--input-file',dest='input_file',default=None,help=help_msg)

        self.args = cmd_parser.parse_args()

    # ----------------------------------------------------------------------------------------------

    def get_input_file(self):
        """
        get input file from cmd line parser
        """
        input_file = self.args.input_file
        
        if input_file is None:
            self.input_file = None
            self.input_path = None
            return

        input_file = os.path.abspath(input_file)
        self.input_path = os.path.dirname(input_file)
        self.input_file = os.path.basename(input_file)

# --------------------------------------------------------------------------------------------------

