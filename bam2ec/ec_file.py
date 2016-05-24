# -*- coding: utf-8 -*-

import logging
import sys
import numpy as np
from collections import OrderedDict

import bam2ec.util as util

LOG = logging.getLogger('BAM2EC')


class EC:
    def __init__(self, filename=None):
        self.version = -1
        self.filename = filename

        self._targets_list = []
        self._targets_dict = OrderedDict()

        self._haplotypes_list = []
        self._haplotypes_dict = OrderedDict()

        # Version 0

        self._reads_list = []
        self._reads_dict = OrderedDict()

        # Version 1

        self._ec_list = []
        self._ec_counts_list = []

        self._alignments = []


def parse(file_in):

    if not file_in:
        raise ValueError("empty file name, cannot load")

    LOG.info("EC File: {0}".format(file_in))

    f = open(file_in, 'rb')

    ec = EC(file_in)

    ec.version = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]

    if ec.version == 0:
        LOG.info("Version: 0, Reads")
    elif ec.version == 1:
        LOG.info("Version: 1, Equivalence Class")
    else:
        LOG.info("Unknown version, exiting")
        LOG.info("Exiting")
        sys.exit()

    # TARGETS

    num_targets = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
    LOG.info("Target Count: {0:,}".format(num_targets))

    ec._targets_list = []
    ec._targets_dict = OrderedDict()

    for i in xrange(0, num_targets):
        str_len = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        target = np.fromfile(f, dtype=np.dtype('a' + str(str_len)), count=1)[0]
        ec._targets_dict[target] = i
        ec._targets_list.append(target)

    LOG.debug('ec._targets_list[0:10]={}'.format(str(ec._targets_list[0:10])))
    LOG.debug('len(ec._targets_list)={}'.format(len(ec._targets_list)))


    # HAPLOTYPES

    num_haplotypes = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
    LOG.info("Haplotype Count: {0:,}".format(num_haplotypes))

    ec._haplotypes_list = []
    ec._haplotypes_dict = OrderedDict()

    for i in xrange(0, num_haplotypes):
        str_len = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        haplotype = np.fromfile(f, dtype=np.dtype('a' + str(str_len)), count=1)[0]
        ec._haplotypes_dict[haplotype] = i
        ec._haplotypes_list.append(haplotype)


    if ec.version == 0:

        # READS

        ec._reads_list = []
        ec._reads_dict = OrderedDict()

        num_reads = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.info("Read Count: {0:,}".format(num_reads))

        for i in xrange(0, num_reads):
            str_len = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
            read_id = np.fromfile(f, dtype=np.dtype('a' + str(str_len)), count=1)[0]
            ec._reads_dict[read_id] = i
            ec._reads_list.append(read_id)

        # ALIGNMENTS

        num_alignments = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.info("Alignment Count: {0:,}".format(num_alignments))

        temp_alignments = np.fromfile(f, dtype = np.dtype('i'), count=num_alignments*3)
        ec._alignments = []

        for i in xrange(0, num_alignments*3, 3):
            read_index = temp_alignments[i]
            target_index = temp_alignments[i+1]
            bit_flag = temp_alignments[i+2]
            ec._alignments.append((read_index, target_index, bit_flag))
    else:


        num_ec = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.info("Equivalance Class Count: {0:,}".format(num_ec))

        ec._ec_list = [x for x in xrange(0, num_ec)]
        ec._ec_counts_list = np.fromfile(f, dtype=np.dtype('i'), count=num_ec)

        # ALIGNMENTS

        num_alignments = np.fromfile(f, dtype=np.dtype('i'), count=1)[0]
        LOG.info("Alignment Count: {0:,}".format(num_alignments))

        temp_alignments = np.fromfile(f, dtype=np.dtype('i'), count=num_alignments*3)
        ec._alignments = []

        for i in xrange(0, num_alignments*3, 3):
            ec_index = temp_alignments[i]
            target_index = temp_alignments[i+1]
            bit_flag = temp_alignments[i+2]
            ec._alignments.append((ec_index, target_index, bit_flag))

    return ec


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

                bits = util.int_to_list(temp_bits, num_haplotypes)
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

                    bits = util.int_to_list(temp_bits, num_haplotypes)
                    LOG.info("{}\t{}\t{}\t".format(rid, target_ids[lid], bits))
    except:
        util._show_error()


