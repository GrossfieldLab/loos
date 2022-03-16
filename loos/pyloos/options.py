#!/usr/bin/env python3

import argparse
import sys


class FullHelper(argparse.Action):
    def __init__(self, option_strings, fullhelp, *args, **kwargs):
        kwargs['nargs'] = 0
        self._fullhelp = fullhelp
        super(FullHelper, self).__init__(option_strings, *args, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        print(self._fullhelp)
        parser.print_help()
        setattr(namespace, self.dest, True)
        parser.exit()


class LoosOptions:
    """
        Parse command line options and implement common sets of options
        Implemented as a wrapper around argparse
    """

    def __init__(self, description, fullhelp=None):
        self.parser = argparse.ArgumentParser(description=description)

        self.optionSets = {}

        if fullhelp:
            self.setFullhelp(fullhelp)

    def setFullhelp(self, fullhelp=None):
        self.parser.add_argument('--fullhelp',
                                 action=FullHelper,
                                 fullhelp=fullhelp)
        self.optionSets['fullhelp'] = True

    # Set up some default arguments
    def modelSelectionPositionalOptions(self):
        self.parser.add_argument(help="Model file describing system contents")
        self.parser.add_argument(help='Use this selection for computation',
                                 default='all')
        self.optionSets['modelSelectionPositionalOptions'] = True

    # Set up some default arguments
    def modelSelectionOptions(self):
        self.parser.add_argument('-m', '--model',
                                 help="Model file describing system contents"
                                 )
        self.parser.add_argument('--selection',
                                 help='Use this selection for computation',
                                 default='all')
        self.optionSets['modelSelectionOptions'] = True

    def trajOptions(self):
        self.parser.add_argument('-t', '--traj',
                                 help='Filename of trajectory or trajectories',
                                 nargs='+')
        self.parser.add_argument('-k', '--skip',
                                 help='Skip frames from the trajectory start',
                                 type=int,
                                 default=0)
        self.parser.add_argument('-s', '--stride',
                                 help='Step through the trajectory by this',
                                 type=int,
                                 default=1)
        self.optionSets['trajOptions'] = True

    def parse_args(self):
        # check this first; otherwise, required arguments will
        # prevent the parser from printing help
        if len(sys.argv) == 1:
            self.parser.print_help(sys.stderr)
            sys.exit(0)

        args = self.parser.parse_args()

        # postprocess args here; error checks, checks for needed behavior
        error = False
        if self.optionSets['modelSelectionOptions']:
            # must supply a model, but since there's a default for selection
            # supplying that is optional
            if not args.model:
                sys.stderr.write('\nRequired option -m/--model missing\n')
                error = True

        if self.optionSets['trajOptions']:
            # must supply a traj, but skip and stride are optional
            if not args.traj:
                sys.stderr.write('\nRequired option -t/--traj missing\n')
                error = True

        if error:
            sys.stderr.write("\n")
            self.parser.print_help(sys.stderr)
            sys.exit(1)

        return args

    def header(self):
        vals = ["'"+str(x) + "'" for x in sys.argv]
        string = " ".join(vals)
        return string
