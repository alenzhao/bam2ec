# -*- coding: utf-8 -*-

import logging
import os
import sys
import traceback

from collections import OrderedDict
from struct import pack

import pysam
import numpy as np

from emase import AlignmentPropertyMatrix as APM

import bam2ec.ec_file as ec_file
import bam2ec.emase_file as emase_file
import bam2ec.util as util

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


"""
--------------------------------------------------------------------
FORMAT                          integer	0 for reads, 1 for equivalence class

NUM_TARGETS	                    integer
    (repeated NUM_TARGETS times)
    TARGET_NAME_LEN 	        integer
    TARGET_NAME	                string, length = TARGET_NAME_LEN

NUM_HAPLOTYPES	                integer
    (repeated NUM_HAPLOTYPES times)
    HAPLOTYPE_NAME_LEN          integer
    HAPLOTYPE_NAME	            string, length = HAPLOTYPE_NAME_LEN

--------------------------------------------------------------------
FORMAT = 0

NUM_READS                       integer
    (repeated NUM_READS times)
    READ_NAME_LEN               integer
    READ_NAME                   string, length = READ_NAME_LEN

NUM_PSEUDO_ALIGN                integer
    (repeated NUM_PSEUDO_ALIGN times)
    READ_INDEX                  integer, index into READS section
    TARGET_INDEX                integer, index into TARGETS section
    BITWISE_FLAG                integer, compressed integer bit field

--------------------------------------------------------------------
FORMAT = 1

NUM_EC                          integer
    (repeated NUM_EC times)
    EC_CLASS_COUNT              integer, the number of times the ec appears

NUM_PSEUDO_ALIGN                integer
    (repeated NUM_PSEUDO_ALIGN times)
    EC_INDEX                    integer, index into EC section
    TARGET_INDEX                integer, index into TARGETS section
    BITWISE_FLAG                integer, compressed integer bit field

--------------------------------------------------------------------

* BITWISE_FLAG is an integer representing the number of haplotypes
that are "turned on" for a target

"""


def convert(file_in, file_out, target_file=None, emase=False):
    """

    :param file_in: Input BAM/SAM file.
    :param file_out: Output file name.
    :param target_file: The target file is a list of main targets that will be used as main targets,
                        not to limit the main targets.  Useful for comparison purposes between BAM files.
    :param emase: Emase output or normal.
    :return:
    """
    LOG.info('Input File: {}'.format(file_in))
    LOG.info('Output File: {}'.format(file_out))

    if target_file:
        LOG.info('Target File: {}'.format(target_file))

    if emase:
        LOG.info('Emase format requested')

    main_targets = OrderedDict()

    if target_file:
        main_targets = util.parse_target_file(target_file)
        if len(main_targets) == 0:
            LOG.error("Unable to parse target file")
            sys.exit(-1)

    # ec = equivalence class
    #      the KEY is a comma separated string of tids
    #      the VALUE is the number of times this equivalence class has appeared
    ec = OrderedDict()

    # ec_idx = lookup to ec
    #          the KEY is a comma separated string of tids
    #          the VALUE is a number specifying the insertion order of the KEY value in ec
    ec_idx = {}

    # all the haplotypes
    haplotypes = set()

    # a lookup of tids to main_targets (Ensembl IDs)
    target_idx_to_main_target = {}

    # unique number of tids encountered and the count
    unique_tids = {}

    # unique reads
    unique_reads = {}

    # times encountering new read id
    read_id_switch_counter = 0

    same_read_target_counter = 0

    try:
        sam_file = pysam.Samfile(file_in, 'rb')
        if len(sam_file.header) == 0:
            raise Exception("BAM File has no header information")
    except:
        sam_file = pysam.Samfile(file_in, 'r')
        if len(sam_file.header) == 0:
            raise Exception("SAM File has no header information")

    line_no = 0
    ec_key = None
    tid = None

    target_ids = []
    try:
        read_id = None

        while True:
            alignment = sam_file.next()
            line_no += 1

            # reference_sequence_name = Column 3 from file, the Reference NAME (EnsemblID_Haplotype)
            # tid = the target id, which is 0 or a positive integer mapping to entries
            #       within the sequence dictionary in the header section of a BAM file
            # main_target = the Ensembl id of the transcript

            reference_sequence_name = sam_file.getrname(alignment.tid)
            tid = str(alignment.tid)
            main_target = reference_sequence_name.split('_')[0]

            try:
                unique_tids[tid] += 1
            except KeyError:
                unique_tids[tid] = 1

            #LOG.verbose("{}\t{}\t{}".format(main_target, reference_sequence_name, tid))

            if target_file:
                if main_target not in main_targets:
                    LOG.error("Unexpected target found in BAM file: {}".format(main_target))
                    sys.exit(-1)
            else:
                if main_target not in main_targets:
                    main_targets[main_target] = len(main_targets)

            target_idx_to_main_target[tid] = main_target
            haplotypes.add(reference_sequence_name.split('_')[1])

            # read_id = Column 1 from file, the Query template NAME
            if read_id is None:
                read_id = alignment.qname

            try:
                unique_reads[read_id] += 1
            except KeyError:
                unique_reads[read_id] = 1

            if read_id != alignment.qname:
                ec_key = ','.join(sorted(target_ids))

                try:
                    ec[ec_key] += 1
                except KeyError:
                    ec[ec_key] = 1
                    ec_idx[ec_key] = len(ec_idx)

                read_id = alignment.qname
                target_ids = [tid]
                read_id_switch_counter += 1
            else:
                if tid not in target_ids:
                    target_ids.append(tid)
                else:
                    same_read_target_counter += 1

            if line_no % 1000000 == 0:
                LOG.info("{0:,} alignments processed, with {1:,} equivalence classes".format(line_no, len(ec)))

    except StopIteration:
        LOG.info("{0:,} alignments processed, with {1:,} equivalence classes".format(line_no, len(ec)))

    if tid not in target_ids:
        target_ids.append(tid)
    else:
        same_read_target_counter += 1

    ec_key = ','.join(sorted(target_ids))

    try:
        ec[ec_key] += 1
    except KeyError:
        ec[ec_key] = 1
        ec_idx[ec_key] = len(ec_idx)

    haplotypes = sorted(list(haplotypes))

    LOG.info("# Unique Reads: {:,}".format(len(unique_reads)))
    LOG.info("# Reads/Target Duplications: {:,}".format(same_read_target_counter))
    LOG.info("# Main Targets: {:,}".format(len(main_targets)))
    LOG.info("# Haplotypes: {:,}".format(len(haplotypes)))
    LOG.info("# Unique Targets: {:,}".format(len(unique_tids)))
    LOG.info("# Equivalence Classes: {:,}".format(len(ec)))

    try:
        os.remove(file_out)
    except OSError:
        pass

    if emase:
        try:
            LOG.info('Creating APM...')
            if LOG.isEnabledFor(util.VERBOSE_LEVELV_NUM):
                LOG.verbose("HAPLOTYPES")
                for h in haplotypes:
                    LOG.verbose(h)
                LOG.verbose("MAIN TARGETS")
                for m in main_targets:
                    LOG.verbose(m)

            new_shape = (len(main_targets), len(haplotypes), len(ec))

            ec_ids = [x for x in xrange(0, len(ec))]

            LOG.debug('Shape={}'.format(new_shape))

            apm = APM(shape=new_shape, haplotype_names=haplotypes, locus_names=main_targets.keys(), read_names=ec_ids)

            # ec.values -> the number of times this equivalence class has appeared
            apm.count = ec.values()

            # k = comma seperated string of tids
            # v = the count
            for k, v in ec.iteritems():
                arr_target_idx = k.split(",")

                # get the main targets by name
                temp_main_targets = set()
                for idx in arr_target_idx:
                    temp_main_targets.add(target_idx_to_main_target[idx])

                # loop through the targets and haplotypes to get the bits
                for main_target in temp_main_targets:
                    # main_target is not an index, but a value like 'ENMUST..001'

                    for i, hap in enumerate(haplotypes):
                        read_transcript = '{}_{}'.format(main_target, hap) # now 'ENMUST..001_A'
                        # get the numerical tid corresponding to read_transcript
                        read_transcript_idx = str(sam_file.gettid(read_transcript))

                        if read_transcript_idx in arr_target_idx:
                            LOG.debug("{}\t{}\t{}".format(ec_idx[k], main_targets[main_target], i))

                            # main_targets[main_target] = idx of main target
                            # i = the haplotype
                            # ec_idx[k] = index of ec
                            apm.set_value(main_targets[main_target], i, ec_idx[k], 1)

            LOG.info("Finalizing...")
            apm.finalize()
            apm.save(file_out, title='bam2ec')
        except:
            util._show_error()
    else:
        try:
            LOG.info("Generating BIN file...")

            f = open(file_out, "wb")

            # version
            f.write(pack('<i', 1))
            LOG.verbose("1\t# VERSION")

            # targets
            LOG.verbose("{:,}\t# NUMBER OF TARGETS".format(len(main_targets)))
            f.write(pack('<i', len(main_targets)))
            for main_target, idx in main_targets.iteritems():
                LOG.verbose("{:,}\t{}\t# {:,}".format(len(main_target), main_target, idx))
                f.write(pack('<i', len(main_target)))
                f.write(pack('<{}s'.format(len(main_target)), main_target))

            # haplotypes
            LOG.verbose("{:,}\t# NUMBER OF HAPLOTYPES".format(len(haplotypes)))
            f.write(pack('<i', len(haplotypes)))
            for idx, hap in enumerate(haplotypes):
                LOG.verbose("{:,}\t{}\t# {:,}".format(len(hap), hap, idx))
                f.write(pack('<i', len(hap)))
                f.write(pack('<{}s'.format(len(hap)), hap))

            # equivalence classes
            LOG.verbose("{:,}\t# NUMBER OF EQUIVALANCE CLASSES".format(len(ec)))
            f.write(pack('<i', len(ec)))
            for idx, k in enumerate(ec.keys()):
                # ec[k] is the count
                LOG.verbose("{:,}\t# {}\t{:,}".format(ec[k], k, idx))
                f.write(pack('<i', ec[k]))

            LOG.info("Determining mappings...")

            # equivalence class mappings
            counter = 0
            for k, v in ec.iteritems():
                arr_target_idx = k.split(",")

                # get the main targets by name
                temp_main_targets = set()
                for idx in arr_target_idx:
                    temp_main_targets.add(target_idx_to_main_target[idx])

                counter += len(temp_main_targets)

            LOG.verbose("{:,}\t# NUMBER OF EQUIVALANCE CLASS MAPPINGS".format(counter))
            f.write(pack('<i', counter))

            for k, v in ec.iteritems():
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

                    LOG.verbose("{}\t{}\t{}\t# {}\t{}".format(ec_idx[k], main_targets[main_target], util.list_to_int(bits), main_target, bits))
                    f.write(pack('<i', ec_idx[k]))
                    f.write(pack('<i', main_targets[main_target]))
                    f.write(pack('<i', util.list_to_int(bits)))

            f.close()
        except:
            util._show_error()

    LOG.info("Done with converting BAM file!")

