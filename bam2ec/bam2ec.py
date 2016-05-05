# -*- coding: utf-8 -*-

import logging
import os
import sys

from struct import pack
from collections import OrderedDict

from emase import AlignmentPropertyMatrix as APM

import pysam

from . import utils


LOG = utils.get_logger()


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
        main_targets = utils.parse_target_file(target_file)
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

    try:
        read_id = None

        target_ids = []
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
                LOG.info("{0:,} lines processed, {1:,} ec classes".format(line_no, len(ec)))

    except StopIteration:
        LOG.info("{0:,} lines processed, {1:,} ec classes".format(line_no, len(ec)))

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
            if LOG.isEnabledFor(utils.VERBOSE_LEVELV_NUM):
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
            utils._show_error()
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

                    LOG.verbose("{}\t{}\t{}\t# {}\t{}".format(ec_idx[k], main_targets[main_target], utils.list_to_int(bits), main_target, bits))
                    f.write(pack('<i', ec_idx[k]))
                    f.write(pack('<i', main_targets[main_target]))
                    f.write(pack('<i', utils.list_to_int(bits)))

            f.close()
        except:
            utils._show_error()

    LOG.info("Done!")


