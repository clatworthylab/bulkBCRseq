#!/usr/bin/env python
from isotyper.utilities._args import COMMAND_MODE  # , REVERSE_PRIMER_GROUP,

# print("Reverse primer group: ", REVERSE_PRIMER_GROUP)

# Commands
if COMMAND_MODE == "1":
    from isotyper.qualitycontrol import preliminary

    preliminary.main()

if COMMAND_MODE == "2":
    from isotyper.qualitycontrol import qualitycontrol

    qualitycontrol.main()


# Clustering reads
if COMMAND_MODE == "3":
    from isotyper.network import network

    network.main()
