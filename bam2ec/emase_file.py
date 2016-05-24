# -*- coding: utf-8 -*-

import logging
from collections import OrderedDict

import emase

LOG = logging.getLogger('BAM2EC')


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


class EMASE:
    def __init__(self, filename=None):
        self.filename = filename

        self._target_list = []
        self._target_dict = OrderedDict()

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

    LOG.info("Emase File: {0}".format(file_in))

    #f = open(file_in, 'rb')

    em = EMASE(file_in)

    apm = emase.AlignmentPropertyMatrix(h5file=file_in)

    em._target_list = list(apm.lname)
    em._target_dict = {target: idx for idx, target in enumerate(em._target_list)}

    em._haplotypes_list = apm.hname
    em._haplotypes_dict = {hap: idx for idx, hap in enumerate(em._haplotypes_list)}

    em._ec_list = list(apm.rname)
    em._ec_counts_list = list(apm.count)

    for ec_idx, ec in enumerate(em._ec_list):
        for target, target_idx in em._target_dict.iteritems():
            bits = []

            for hap, hap_idx in em._haplotypes_dict.iteritems():
                bits.append(apm.data[hap_idx][ec_idx, target_idx])

            em._alignments.append((ec_idx, target_idx, list_to_int(bits)))

    return em






