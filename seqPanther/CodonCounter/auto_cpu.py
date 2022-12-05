#!/usr/bin/env python
"""Assigns number of CPU should be use by the application."""

from multiprocessing import cpu_count, Pool
import warnings


def cpus(ncr):
    """Return number of cpus to be used based on user provided cores.
    uses all the core is 0 is provided.
    if given number of cores is greater than the number of cores available,
    then it will use n-2 cores. where n is the number of cores available.

    :ncr: number of cores to be used.
    :returns: pool

    """
    if ncr <= 0 or ncr > cpu_count():
        ncr = cpu_count() - 2
        warnings.warn("Number of cores provided is either 0, negative or more "
                      "than the number of cores available.\n"
                      "Using {} cores".format(ncr))
    return Pool(ncr)
