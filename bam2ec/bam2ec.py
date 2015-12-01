# -*- coding: utf-8 -*-

import logging
import traceback
from struct import pack
from collections import OrderedDict
import sys

import sys
import time
import argparse
import numpy as np
from collections import OrderedDict
from emase import AlignmentPropertyMatrix as APM

import pysam


LOG = logging.getLogger('BAM2EC')

#
#
# Logging
#
#


class BAM2ECFormatter(logging.Formatter):
    err_fmt = "[bam2ec] %(msg)s"
    # err_fmt = "DEBUG: %(asctime)s, %(module)s: %(lineno)d: %(msg)s"
    # dbg_fmt = "DEBUG: %(asctime)s, %(module)s: %(lineno)d: %(msg)s"
    dbg_fmt = "[bam2ec debug] %(msg)s"
    info_fmt = "[bam2ec] %(msg)s"

    def __init__(self, fmt="%(levelno)s: %(msg)s"):
        logging.Formatter.__init__(self, fmt)

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._fmt = BAM2ECFormatter.dbg_fmt

        elif record.levelno == logging.INFO:
            self._fmt = BAM2ECFormatter.info_fmt

        elif record.levelno == logging.ERROR:
            self._fmt = BAM2ECFormatter.err_fmt

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._fmt = format_orig

        return result


def get_logger():
    """
    Get the logger
    :return: logger object
    """
    return LOG


def configure_logging(level):
    """

    :param level: 0=ERROR, 1=INFO, 2+=DEBUG
    """
    global LOG
    LOG = logging.getLogger('BAM2EC')

    handler = logging.StreamHandler(sys.stderr)

    mfmt = BAM2ECFormatter()
    handler.setFormatter(mfmt)

    if level == 0:
        LOG.setLevel(logging.INFO)
    else:
        LOG.setLevel(logging.DEBUG)

    LOG.addHandler(handler)


def int_to_list(c, size):
    ret = [0]*size
    for i in xrange(0, size):
        if (c & (1 << i)) != 0:
            ret[i] = 1
    return ret


def list_to_int(l):
    c = 0
    for i, onoff in enumerate(l):
        if onoff == 1:
            c |= 1 << i
    return c

def _show_error():
    """
    show system errors
    """
    et, ev, tb = sys.exc_info()

    print "Error Type: %s" % et
    print "Error Value: %s" % ev
    while tb :
        co = tb.tb_frame.f_code
        filename = str(co.co_filename)
        line_no = str(traceback.tb_lineno(tb))
        print '    %s:%s' % (filename, line_no)
        tb = tb.tb_next


def parse_target_file(target_file):
    targets = OrderedDict()
    with open(target_file, 'r') as f:
        for line in f:
            if line and line[0] == '#':
                continue
            _id = line.strip().split()[0]
            targets[_id] = len(targets)
    return targets


def convert(file_in, file_out, target_file=None, emase=False):
    LOG.info('Input File: {}'.format(file_in))
    LOG.info('Output File: {}'.format(file_out))

    if target_file:
        LOG.info('Target File: {}'.format(target_file))

    if emase:
        LOG.info('Emase format requested')

    main_targets = OrderedDict()

    if target_file:
        main_targets = parse_target_file(target_file)
        if len(main_targets) == 0:
            LOG.error("Unable to parse target file")
            sys.exit(-1)

    ec = OrderedDict()
    ec_idx = {}

    haplotypes = set()

    target_idx_to_main_target = {}

    try:
        sam_file = pysam.Samfile(file_in, 'rb')
        if len(sam_file.header) == 0:
            raise Exception("BAM File has no header information")
    except:
        sam_file = pysam.Samfile(file_in, 'r')
        if len(sam_file.header) == 0:
            raise Exception("SAM File has no header information")

    try:
        line_no = 0
        read_id = None

        temp_transcripts = []
        while True:
            line_no += 1
            alignment = sam_file.next()
            read_transcript = sam_file.getrname(alignment.tid)
            read_transcript_idx = str(alignment.tid)
            #print read_transcript, read_transcript_idx

            main_target = read_transcript.split('_')[0]

            if target_file:
                if main_target not in main_targets:
                    LOG.error("Unexpected target found in BAM file: {}".format(main_target))
                    sys.exit(-1)

            else:
                if main_target not in main_targets:
                    main_targets[main_target] = len(main_targets)

            target_idx_to_main_target[read_transcript_idx] = main_target
            haplotypes.add(read_transcript.split('_')[1])

            if read_id is None:
                read_id = alignment.qname

            if read_id != alignment.qname:
                ec_key = ','.join(sorted(temp_transcripts))

                try:
                    ec[ec_key] += 1
                except KeyError:
                    ec[ec_key] = 1
                    ec_idx[ec_key] = len(ec_idx)

                read_id = alignment.qname
                temp_transcripts = [read_transcript_idx]
            else:
                temp_transcripts.append(read_transcript_idx)

            if line_no % 1000000 == 0:
                LOG.info("{0:,} reads processed, with {1:,} ec classes".format(line_no, len(ec)))

    except StopIteration:
        LOG.info("{0:,} reads processed, with {1:,} ec classes".format(line_no, len(ec)))

    haplotypes = sorted(list(haplotypes))

    if emase:
        try:
            LOG.info('Creating APM...')
            new_shape = (len(main_targets), len(haplotypes), len(ec))

            ec_ids = [x for x in xrange(0, len(ec))]

            LOG.debug('Shape={}'.format(new_shape))

            aln_mat_kallisto = APM(shape=new_shape, haplotype_names=haplotypes, locus_names=main_targets.keys(), read_names=ec_ids)

            aln_mat_kallisto.count = ec.values()

            for k,v in ec.iteritems():
                arr_target_idx = k.split(",")

                # get the main targets by name
                temp_main_targets = set()
                for idx in arr_target_idx:
                    temp_main_targets.add(target_idx_to_main_target[idx])

                # loop through the haplotypes and targets to get the bits
                for main_target in temp_main_targets:
                    # main_target is not an index, but a value like 'ENMUST..001'

                    for i, hap in enumerate(haplotypes):
                        read_transcript = '{}_{}'.format(main_target, hap) # now 'ENMUST..001_A'
                        read_transcript_idx = str(sam_file.gettid(read_transcript))

                        if read_transcript_idx in arr_target_idx:
                            LOG.debug("{}\t{}\t{}".format(ec_idx[k], main_targets[main_target], i))
                            aln_mat_kallisto.set_value(main_targets[main_target], i, ec_idx[k], 1)


            LOG.info("Finalizing...")
            aln_mat_kallisto.finalize()
            aln_mat_kallisto.save(file_out, title='bam2ec')

            LOG.info("DONE")
        except:
            _show_error()


    else:
        LOG.info("Generating file...")
        f = open("{}.bin".format(file_in), "wb")

        # version
        f.write(pack('<i', 1))

        # targets
        f.write(pack('<i', len(main_targets)))
        for main_target, idx in main_targets.iteritems():
            f.write(pack('<i', len(main_target)))
            f.write(pack('<{}s'.format(len(main_target)), main_target))

        # haplotypes
        f.write(pack('<i', len(haplotypes)))
        for hap in haplotypes:
            f.write(pack('<i', len(hap)))
            f.write(pack('<{}s'.format(len(hap)), hap))

        # ec
        f.write(pack('<i', len(ec)))
        for k,v in ec.iteritems():
            # v is the count
            f.write(pack('<i', v))

        LOG.info("Determining mappings...")

        # ec mappings
        counter = 0
        for k,v in ec.iteritems():
            arr_target_idx = k.split(",")

            # get the main targets by name
            temp_main_targets = set()
            for idx in arr_target_idx:
                temp_main_targets.add(target_idx_to_main_target[idx])

            counter += len(temp_main_targets)

        f.write(pack('<i', counter))
        for k,v in ec.iteritems():
            arr_target_idx = k.split(",")

            # get the main targets by name
            temp_main_targets = set()
            for idx in arr_target_idx:
                temp_main_targets.add(target_idx_to_main_target[idx])

            # loop through the haplotypes and targets to get the bits
            for main_target in temp_main_targets:
                # main_target is not an index, but a value like 'ENMUST..001'

                bits = []

                for hap in haplotypes:
                    read_transcript = '{}_{}'.format(main_target, hap) # now 'ENMUST..001_A'
                    read_transcript_idx = str(sam_file.gettid(read_transcript))

                    if read_transcript_idx in arr_target_idx:
                        bits.append(1)
                    else:
                        bits.append(0)

                LOG.debug("{}\t{}\t{}\t{}".format(ec_idx[k], main_target, bits, list_to_int(bits)))
                f.write(pack('<i', ec_idx[k]))
                f.write(pack('<i', main_targets[main_target]))
                f.write(pack('<i', list_to_int(bits)))

        f.close()

    LOG.info("Done!")


def view(binary_file_name, detail=False):
    """

    :param binary_file_name:
    :return:
    """

    if not binary_file_name:
        raise ValueError("empty file name, cannot load")

    LOG.info("Binary File: {0}".format(binary_file_name))

    f = open(binary_file_name, 'rb')

    file_version = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]

    if file_version == 0:
        LOG.info("Version: 0, Reads")
    elif file_version == 1:
        LOG.info("Version: 1, Equivalence Class")
    else:
        LOG.info("Unknown version, exiting")

    # TARGETS

    target_ids = []
    targets = OrderedDict()

    num_targets = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
    LOG.info("Target Count: {0:,}".format(num_targets))

    for i in xrange(0, num_targets):
        str_len = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        target = np.fromfile(f, dtype=np.dtype('a' + str(str_len)), count=1)[0]
        targets[target] = i
        target_ids.append(target)
        if detail:
            LOG.info("{}\t{}".format(i, target))

    # HAPLOTYPES

    haplotype_ids = []
    haplotypes = OrderedDict()

    num_haplotypes = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
    LOG.info("Haplotype Count: {0:,}".format(num_haplotypes))

    for i in xrange(0, num_haplotypes):
        str_len = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        haplotype = np.fromfile(f, dtype=np.dtype('a' + str(str_len)), count=1)[0]
        haplotypes[haplotype] = i
        haplotype_ids.append(haplotype)
        if detail:
            LOG.info("{}\t{}".format(i, haplotype))

    if file_version == 0:

        # READS

        read_ids = []
        reads = OrderedDict()

        num_reads = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.info("Read Count: {0:,}".format(num_reads))

        for i in xrange(0, num_reads):
            str_len = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
            read_id = np.fromfile(f, dtype=np.dtype('a' + str(str_len)), count=1)[0]
            reads[read_id] = i
            read_ids.append(read_id)

        # ALIGNMENTS

        num_alignments = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.info("Alignment Count: {0:,}".format(num_alignments))

        alignments = np.fromfile(f, dtype = np.dtype('i'), count=num_alignments*3)

        counter = 0

        for i in xrange(0, num_alignments*3, 3):
            rid = alignments[i]
            lid = alignments[i+1]
            temp_bits = alignments[i+2]

            counter += 1
            if temp_bits == 0:
                continue

            bits = int_to_list(temp_bits, num_haplotypes)
            if detail:
                print rid, target_ids[lid], bits

    else:

        # EQUIVALENCE CLASSES

        num_ec = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.info("Equivalance Class Count: {0:,}".format(num_ec))

        ec_ids = [x for x in xrange(0, num_ec)]
        counts = np.fromfile(f, dtype=np.dtype('i'), count=num_ec)

        # ALIGNMENTS

        num_alignments = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.info("Alignment Count: {0:,}".format(num_alignments))

        if detail:
            alignments = np.fromfile(f, dtype=np.dtype('i'), count=num_alignments*3)

            counter = 0

            for i in xrange(0, num_alignments*3, 3):
                rid = alignments[i]
                lid = alignments[i+1]
                temp_bits = alignments[i+2]

                counter += 1
                if temp_bits == 0:
                    continue

                bits = int_to_list(temp_bits, num_haplotypes)
                LOG.info("{}\t{}\t{}\t".format(rid, target_ids[lid], bits))
