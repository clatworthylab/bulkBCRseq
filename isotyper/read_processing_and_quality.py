#!/usr/bin/env python
from isotyper.utilities._args import COMMAND_MODE

# prepare for QC
if COMMAND_MODE == "1":
    from isotyper.qualitycontrol import preliminary

    preliminary.main()

# Main QC
if COMMAND_MODE == "2":
    from isotyper.qualitycontrol import qualitycontrol

    qualitycontrol.main()

# Clustering reads
if COMMAND_MODE == "3":
    from isotyper.network import network

    network.main()

# Clustering reads
if COMMAND_MODE == "4":
    from isotyper.statistics import statistics

    statistics.main()
