# -*- coding: utf-8 -*-

import logging
import sys
import traceback

from collections import OrderedDict
from struct import pack

import numpy as np
from emase import AlignmentPropertyMatrix as APM

from . import ec_file
from . import emase_file

VERBOSE_LEVELV_NUM = 9


def verbose(self, message, *args, **kws):
    # Yes, logger takes its '*args' as 'args'.
    if self.isEnabledFor(VERBOSE_LEVELV_NUM):
        self._log(VERBOSE_LEVELV_NUM, message, args, **kws)


logging.addLevelName(VERBOSE_LEVELV_NUM, "DEBUGV")
logging.Logger.verbose = verbose

LOG = logging.getLogger('BAM2EC')



#
#
# Logging
#
#
# CRITICAL  50
# ERROR     40
# WARNING   30
# INFO      20
# DEBUG     10  -d
# VERBOSE   9   -dd
# NOTSET    0
#
#


class BAM2ECFormatter(logging.Formatter):
    dbg_fmt = "[bam2ec][%(asctime)s] %(msg)s"
    info_fmt = "[bam2ec] %(msg)s"
    # err_fmt = "DEBUG: %(asctime)s, %(module)s: %(lineno)d: %(msg)s"
    # dbg_fmt = "DEBUG: %(asctime)s, %(module)s: %(lineno)d: %(msg)s"

    def __init__(self, fmt="%(asctime)s, %(module)s: %(lineno)d: %(msg)s"):
        logging.Formatter.__init__(self, fmt)

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._fmt

        # Replace the original format with one customized by logging level
        if record.levelno <= logging.DEBUG:
            self._fmt = BAM2ECFormatter.dbg_fmt

        elif record.levelno == logging.INFO:
            self._fmt = BAM2ECFormatter.info_fmt

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
    elif level == 1:
        LOG.setLevel(logging.DEBUG)
    elif level >= 2:
        LOG.setLevel(VERBOSE_LEVELV_NUM)

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


def dump(binary_file_name, detail=False):
    """

    :param binary_file_name:
    :return:
    """

    try:

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
    except:
        _show_error()



def bin2emase(binary_file_name, emase_file_name):
    try:
        if not binary_file_name:
            raise ValueError("empty file name, cannot load")

        LOG.info("Binary File: {0}".format(binary_file_name))

        f = open(binary_file_name, 'rb')

        file_version = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]

        if file_version == 0:
            LOG.info("Version: 0, Reads, exiting")
            sys.exit(-1)
        elif file_version == 1:
            LOG.info("Version: 1, Equivalence Class")
        else:
            LOG.info("Unknown version, exiting")
            sys.exit(-1)

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

            LOG.verbose("{}\t{}".format(i, target))

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

            LOG.verbose("{}\t{}".format(i, haplotype))

        # EQUIVALENCE CLASSES

        num_ec = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.info("Equivalance Class Count: {0:,}".format(num_ec))

        ec_ids = [x for x in xrange(0, num_ec)]
        counts = np.fromfile(f, dtype=np.dtype('i'), count=num_ec)

        # ALIGNMENTS

        num_alignments = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.info("Alignment Count: {0:,}".format(num_alignments))


        new_shape = (num_targets, num_haplotypes, num_ec)
        LOG.info('Creating APM...')

        ec_ids = [x for x in xrange(0, num_ec)]

        LOG.debug('Shape={}'.format(new_shape))

        apm = APM(shape=new_shape, haplotype_names=haplotype_ids, locus_names=target_ids, read_names=ec_ids)

        # counts -> the number of times this equivalence class has appeared
        apm.count = counts

        counter = 0
        alignments = np.fromfile(f, dtype=np.dtype('i'), count=num_alignments*3)

        for i in xrange(0, num_alignments*3, 3):
            rid = alignments[i]
            lid = alignments[i+1]
            temp_bits = alignments[i+2]

            counter += 1
            if temp_bits == 0:
                continue

            try:
                bits = int_to_list(temp_bits, num_haplotypes)
                for i, bit in enumerate(bits):
                    if bit:
                        apm.set_value(rid, i, lid, 1)
            except Exception, e:
                print rid
                print i
                print lid
                _show_error()
                raise e

        LOG.info("Finalizing...")
        apm.finalize()
        apm.save(emase_file_name, title='bam2ec')

    except:
        _show_error()


def ec2emase(file_in, file_out, target_file=None):
    ec = ec_file.parse(file_in)
    new_shape = (len(ec._targets_list), len(ec._haplotypes_list), len(ec._ec_list))

    LOG.info('Creating APM...')
    LOG.debug('Shape={}'.format(new_shape))


    apm = APM(shape=new_shape, haplotype_names=ec._haplotypes_list, locus_names=ec._targets_list, read_names=ec._ec_list)

    LOG.debug('ec._haplotypes_list={}'.format(str(ec._haplotypes_list)))
    LOG.debug('ec._targets_list[0:10]={}'.format(str(ec._targets_list[0:10])))
    LOG.debug('ec._ec_list[0:10]={}'.format(str(ec._ec_list[0:10])))

    # counts -> the number of times this equivalence class has appeared
    apm.count = ec._ec_counts_list

    counter = 0
    num_haplotypes = len(ec._haplotypes_list)

    try:

        for alignment in ec._alignments:
            #LOG.verbose(str(alignment))
            ec_index = alignment[0]
            target_index = alignment[1]
            temp_bits = alignment[2]

            if temp_bits == 0:
                continue

            bits = int_to_list(temp_bits, num_haplotypes)
            for i, bit in enumerate(bits):
                if bit:
                    # lid, hid, rid, value
                    apm.set_value(target_index, i, ec_index, 1)
    except Exception, e:
        _show_error()
        raise e

    LOG.info("Finalizing...")
    apm.finalize()
    apm.save(file_out, title='bam2ec')





def emase2ec(file_in, file_out):
    emase = emase_file.parse(file_in)

    try:
        LOG.info("Generating BIN file...")

        f = open(file_out, "wb")

        # version
        f.write(pack('<i', 1))
        LOG.verbose("1\t# VERSION")

        # targets
        LOG.verbose("{:,}\t# NUMBER OF TARGETS".format(len(emase._target_list)))
        f.write(pack('<i', len(emase._target_list)))
        for main_target, idx in emase._target_dict.iteritems():
            LOG.verbose("{:,}\t{}\t# {:,}".format(len(main_target), main_target, idx))
            f.write(pack('<i', len(main_target)))
            f.write(pack('<{}s'.format(len(main_target)), main_target))

        # haplotypes
        LOG.verbose("{:,}\t# NUMBER OF HAPLOTYPES".format(len(emase._haplotypes_list)))
        f.write(pack('<i', len(emase._haplotypes_list)))
        for idx, hap in enumerate(emase._haplotypes_list):
            LOG.verbose("{:,}\t{}\t# {:,}".format(len(hap), hap, idx))
            f.write(pack('<i', len(hap)))
            f.write(pack('<{}s'.format(len(hap)), hap))

        # equivalence classes
        LOG.verbose("{:,}\t# NUMBER OF EQUIVALANCE CLASSES".format(len(emase._ec_list)))
        f.write(pack('<i', len(emase._ec_list)))
        for idx, k in enumerate(emase._ec_counts_list):
            # k is the count
            LOG.verbose("{:,}\t# {:,}".format(k, idx))
            f.write(pack('<i', k))

        LOG.info("Determining mappings...")

        # equivalence class mappings
        LOG.verbose("{:,}\t# NUMBER OF EQUIVALANCE CLASS MAPPINGS".format(len(emase._alignments)))
        f.write(pack('<i', len(emase._alignments)))

        for idx, alignment in enumerate(emase._alignments):
            LOG.verbose("{}\t{}\t{}\t# {}\t{}".format(alignment[0], alignment[1], list_to_int(alignment[2]), emase._target_list[alignment[1]], alignment[2]))
            f.write(pack('<i', alignment[0]))
            f.write(pack('<i', alignment[1]))
            f.write(pack('<i', list_to_int(alignment[2])))

        f.close()
    except:
        _show_error()


