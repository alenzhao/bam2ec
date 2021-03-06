#!/usr/bin/env python

import argparse
import importlib
import os
import sys

from . import __version__ as version
from . import commands

ext_modules = ['pysam', 'emase', 'numpy']
failed_modules = []

logo_text = """

@@@@@@@    @@@@@@   @@@@@@@@@@    @@@@@@   @@@@@@@@   @@@@@@@
@@@@@@@@  @@@@@@@@  @@@@@@@@@@@  @@@@@@@@  @@@@@@@@  @@@@@@@@
@@!  @@@  @@!  @@@  @@! @@! @@!       @@@  @@!       !@@
!@   @!@  !@!  @!@  !@! !@! !@!      @!@   !@!       !@!
@!@!@!@   @!@!@!@!  @!! !!@ @!@     !!@    @!!!:!    !@!
!!!@!!!!  !!!@!!!!  !@!   ! !@!    !!:     !!!!!:    !!!       v""" + version + """
!!:  !!!  !!:  !!!  !!:     !!:   !:!      !!:       :!!
:!:  !:!  :!:  !:!  :!:     :!:  :!:       :!:       :!:
 :: ::::  ::   :::  :::     ::   :: :::::   :: ::::   ::: :::
:: : ::    :   : :   :      :    :: : :::  : :: ::    :: :: :

"""

for dependency in ext_modules:
    try:
        importlib.import_module(dependency)
    except ImportError, ie:
        failed_modules.append(dependency)

if len(failed_modules) > 0:
    sys.stderr.write('Error: The following modules need to be installed: ')
    sys.stderr.write('\t' + ', '.join(failed_modules))
    sys.exit(1)


class BAM2ECToolsApp(object):
    """
    The most commonly used commands are:
       convert       convert file
       dump          view file
       ec2emase      convert binary file to EMASE format
       emase2ec      convert EMASE format to binary file

    """

    def __init__(self):
        self.script_name = os.path.basename(__file__)
        parser = argparse.ArgumentParser(add_help=False)

        def print_message(message):
            sys.stderr.write(logo_text)
            sys.stderr.write('\n')
            sys.stderr.write(message)
            sys.stderr.write('\n')
            sys.stderr.write(BAM2ECToolsApp.__doc__)
            sys.stderr.write('\n')
            sys.exit(1)

        parser.error = print_message

        parser.add_argument('command', nargs='?', help='Subcommand to run')
        parser.add_argument("-h", "--help", dest="help", action="store_true")
        parser.add_argument("-v", "--version", dest="version", action="store_true")

        if len(sys.argv) == 1:
            parser.print_help()
            sys.exit()

        # parse_args defaults to [1:] for args, but need to exclude
        # the rest of the args too, or validation will fail
        args = parser.parse_args(sys.argv[1:2])

        if args.version:
            print version
            sys.exit(1)

        if args.help:
            parser.print_help()
            sys.exit()

        if not args.command:
            parser.print_help()
            sys.exit()

        if not hasattr(self, args.command):
            parser.print_help()
            sys.exit()

        # use dispatch pattern to invoke method with same name
        getattr(self, args.command)()

    def convert(self):
        commands.command_convert(sys.argv[2:], self.script_name + ' convert')

    def dump(self):
        commands.command_dump(sys.argv[2:], self.script_name + ' dump')

    def ec2emase(self):
        commands.command_ec2emase(sys.argv[2:], self.script_name + ' ec2emase')

    def emase2ec(self):
        commands.command_emase2ec(sys.argv[2:], self.script_name + ' emase2ec')

    def logo(self):
        print logo_text

if __name__ == '__main__':
    BAM2ECToolsApp()

