# -*- coding: utf-8 -*-

import argparse
import sys
import time

from . import bam2ec

LOG = bam2ec.get_logger()


def command_convert(raw_args, prog=None):
    """
    Convert a BAM/SAM file

    Usage: convert [-options] -o <Output file> -i <BAM file>

    Required Parameters:
    -i, --input <BAM file>           input file to convert
    -o, --output <output file>       file to create

    Optional Parameters:
    -e, --emase                      Emase file format
    -t, --target <Target file>       target file name

    Help Parameters:
    -h, --help                       print the help and exit
    -d, --debug                      turn debugging on, list multiple times for more messages

    """

    if prog:
        parser = argparse.ArgumentParser(prog=prog, add_help=False)
    else:
        parser = argparse.ArgumentParser(add_help=False)

    def print_message(message=None):
        if message:
            sys.stderr.write(message)
        else:
            sys.stderr.write(command_convert.__doc__)
        sys.stderr.write('\n')
        sys.exit(1)

    parser.error = print_message

    # required
    parser.add_argument("-i", "--input", dest="input", metavar="Input_File")
    parser.add_argument("-o", "--output", dest="output", metavar="Output_File")

    # optional
    parser.add_argument("-e", "--emase", dest="emase", action='store_true')
    parser.add_argument("-t", "--target", dest="target", metavar="Target_File")


    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    bam2ec.configure_logging(args.debug)

    if args.help:
        print_message()

    if not args.input:
        LOG.error("No input file was specified.")
        print_message()

    if not args.output:
        LOG.error("No output file was specified.")
        print_message()

    try:
        bam2ec.convert(args.input, args.output, args.target, args.emase)
    except KeyboardInterrupt, ki:
        LOG.debug(ki)
    except Exception, e:
        LOG.error(e)


def command_view(raw_args, prog=None):
    """
    View the binary formatted file

    Usage: view [-options] -i <BINary file>

    Required Parameters:
    -i, --input <BINary file>        input file

    Optional Parameters:
    -v, --verbose                    verbose output

    Help Parameters:
    -h, --help                       print the help and exit
    -d, --debug                      turn debugging on, list multiple times for more messages

    """

    if prog:
        parser = argparse.ArgumentParser(prog=prog, add_help=False)
    else:
        parser = argparse.ArgumentParser(add_help=False)

    def print_message(message=None):
        if message:
            sys.stderr.write(message)
        else:
            sys.stderr.write(command_view.__doc__)
        sys.stderr.write('\n')
        sys.exit(1)

    parser.error = print_message

    # required
    parser.add_argument("-i", "--input", dest="input", metavar="Input_File")

    # optional
    parser.add_argument("-v", "--verbose", dest="verbose", action='store_true')

    # debugging and help
    parser.add_argument("-h", "--help", dest="help", action='store_true')
    parser.add_argument("-d", "--debug", dest="debug", action="count", default=0)

    args = parser.parse_args(raw_args)

    bam2ec.configure_logging(args.debug)

    if args.help:
        print_message()

    if not args.input:
        LOG.error("No input file was specified.")
        print_message()

    try:
        bam2ec.view(args.input, args.verbose)
    except KeyboardInterrupt, ki:
        LOG.debug(ki)
    except Exception, e:
        LOG.error(e)
