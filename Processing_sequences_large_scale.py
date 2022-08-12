#!/usr/bin/env python
import sys
import os
from pathlib import Path

# import shutil

Path("logs").mkdir(exist_ok=True)


def Get_info(file):
    """Summary

    Parameters
    ----------
    file : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    (
        id,
        sample,
        info,
        gene,
        directory,
        pair,
        pair_final,
        dir,
        platform,
        spec,
        primer,
        other,
        reverse_primer_group,
    ) = ([], [], [], [], [], [], [], [], [], [], [], [], [])
    fh = open(file, "r")
    ind = 0
    for l in fh:
        ind = ind + 1
        if l[0] != "#" and len(l) > 3:
            l = l.replace(
                "/lustre/scratch108/viruses/rbr1/",
                "/lustre/scratch118/infgen/team146/rbr1/",
            )
            l = l.strip().split()
            id.append(l[0])
            sample.append(l[1])
            info.append(l[2])
            gene.append(l[3])
            directory.append(l[4])
            pair.append(l[5])
            pair_final.append(l[6])
            dir.append(l[7])
            platform.append(l[8])
            spec.append(l[9])
            if len(l) >= 11:
                primer.append(l[10])
            else:
                primer.append("LIBRARY/FR1_primers.txt")
            if len(l) >= 14:
                reverse_primer_group.append(l[13])
            else:
                reverse_primer_group.append("STANDARD")
            if len(l) >= 12:
                lis = []
                for i in range(11, len(l)):
                    lis.append(l[i])
                other.append(",".join(lis))
            else:
                other.append("")
    fh.close()
    return (
        id,
        sample,
        info,
        gene,
        directory,
        pair,
        pair_final,
        dir,
        platform,
        spec,
        primer,
        other,
        reverse_primer_group,
    )


def Get_info_concat(file):
    """Summary

    Parameters
    ----------
    file : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    (id, sample, directory, downsample) = ([], [], [], [])
    fh = open(file, "r")
    # ind = 0
    for l in fh:
        lx = l.split("\t")
        if len(lx) > 3:
            id.append(lx[0])
            sample.append([lx[1]])
            directory.append(lx[2].strip("\n"))
            downsample.append(lx[3])
        elif len(lx) > 1:
            id.append(lx[0])
            sample.append([lx[1]])
            directory.append(lx[2].strip("\n"))
        else:
            break
    fh.close()
    if len(downsample) > 0:
        return (id, sample, directory, downsample)
    else:
        return (id, sample, directory)


args = sys.argv
# mem = '-R"select[mem>1800] rusage[mem=1800]" -M1800 '
# mem = '-R"select[mem>3800] rusage[mem=3800]" -M3800 '
# mem = '-R"select[mem>8000] rusage[mem=8000]" -M8000 '
mem = '-R"select[mem>16000] rusage[mem=16000]" -M16000 '
# mem = ''
# queue = '-q yesterday'
queue = "-q normal"
# queue = "-q parallel"
# queue= "-q long"
# queue ="-q small"
# queue = "-q test"
# command = 'bsub -G team146 '+queue+' python Processing_sequences.prog '
# span = '-n10 -R"span[hosts=1]"'
span = ""

if len(args) < 5:
    queue = "-q normal"
    command = (
        "python Processing_sequences_large_scale.py [sample file list] [commands (comma separated list)] ",
        "[bsub command: Y/N] [print commands: Y/N] [run commands: Y/N]",
    )
    print("SEQUENCE ANALYSIS PIPELINE: Creates networks from MiSeq data")
    print("USAGE:")
    print(command, "\n")
    os.system("cat Command_outline.txt")
    print("\n")
elif len(args) == 7:
    file = args[1]
    concat_file = args[2]
    command = args[3]
    bsub_command = args[4]
    print_command = args[5]
    run_command = args[6]
    command = command.split(",")
    if bsub_command not in ["Y", "N"]:
        print(
            (
                "python Processing_sequences.py [sample file list] [concat file list] [commands (comma separated list)] ",
                "[bsub command: Y/N] [print commands: Y/N] [run commands: Y/N]\n\tError: bsub command must be: Y or N",
            )
        )
    if bsub_command not in ["Y", "N"]:
        print(
            (
                "python Processing_sequences.py [sample file list] [concat file list] [commands (comma separated list)] ",
                "[bsub command: Y/N] [print commands: Y/N] [run commands: Y/N] \n\tError: print command must be: Y or N",
            )
        )
    if bsub_command not in ["Y", "N"]:
        print(
            (
                "python Processing_sequences.py [sample file list] [concat file list] [commands (comma separated list)] ",
                "[bsub command: Y/N] [print commands: Y/N] [run commands: Y/N] \n\tError: run command must be: Y or N",
            )
        )
    try:
        (ids_, samples_, directories_, downsample_) = Get_info_concat(
            concat_file
        )
    except:
        (ids_, samples_, directories_) = Get_info_concat(concat_file)
    commands = []
    for i in range(0, len(samples_)):
        if "downsample_" in locals():
            id_, sample_, directory_, dwnsamp_ = (
                ids_[i],
                samples_[i],
                directories_[i],
                downsample_[i],
            )
        else:
            id_, sample_, directory_ = ids_[i], samples_[i], directories_[i]
        bsub = ""
        if "3.5" in command:
            span = '-n24 -R"span[hosts=1]"'
        elif "3.51" in command:
            span = '-n8 -R"span[hosts=1]"'
        else:
            span = ""
        if bsub_command == "Y":
            bsub = (
                "bsub -P team205 -G teichlab "
                + queue
                + " -o logs/out_MERGING_"
                + id_
                + " -J "
                + id_
                + " "
                + mem
                + " "
                + span
                + " "
            )
        if "3.5" in command:
            command1 = (
                "python "
                + "isotyper/Network_generation_from_fully_reduced_fasta.py "
                + directory_
                + " "
                + id_
                + " "
                + " ".join(sample_)
            )
            commands.append(bsub + command1)
        if "3.51" in command:
            if "dwnsamp_" in locals():
                command1 = (
                    "python "
                    + "isotyper/Network_generation_from_fully_reduced_fasta_V2.py "
                    + directory_
                    + " "
                    + id_
                    + " "
                    + " ".join(sample_)
                    + " "
                    + str(dwnsamp_)
                )
            else:
                command1 = (
                    "python "
                    + "isotyper/Network_generation_from_fully_reduced_fasta_V2.py "
                    + directory_
                    + " "
                    + id_
                    + " "
                    + " ".join(sample_)
                )
            commands.append(bsub + command1)
    for comm in commands:
        if print_command == "Y":
            print(comm, "\n")
        if run_command == "Y":
            os.system(comm)
else:
    file = args[1]
    command = args[2]
    bsub_command = args[3]
    print_command = args[4]
    run_command = args[5]
    command = command.split(",")
    if bsub_command not in ["Y", "N"]:
        print(
            (
                "python Processing_sequences.py [sample file list] [commands (comma separated list)] ",
                "[bsub command: Y/N] [print commands: Y/N] [run commands: Y/N]\n\tError: bsub command must be: Y or N",
            )
        )
    if bsub_command not in ["Y", "N"]:
        print(
            (
                "python Processing_sequences.py [sample file list] [commands (comma separated list)] ",
                "[bsub command: Y/N] [print commands: Y/N] [run commands: Y/N] \n\tError: print command must be: Y or N",
            )
        )
    if bsub_command not in ["Y", "N"]:
        print(
            (
                "python Processing_sequences.py [sample file list] [commands (comma separated list)] ",
                "[bsub command: Y/N] [print commands: Y/N] [run commands: Y/N] \n\tError: run command must be: Y or N",
            )
        )
    (
        ids,
        samples,
        infos,
        gene,
        source,
        pairs,
        pairs_final,
        dirs,
        platform,
        spec,
        primer,
        others,
        reverse_primer_group,
    ) = Get_info(file)
    print(dir)
    commands = []
    for i in range(0, len(samples)):
        (
            info,
            sample,
            gene_types,
            pair,
            pair_final,
            dir,
            platforms,
            primers,
            other,
        ) = (
            infos[i],
            samples[i],
            gene[i],
            pairs[i],
            pairs_final[i],
            dirs[i],
            platform[i],
            primer[i],
            others[i],
        )
        bsub = ""
        if "3" in command:
            span = '-n10 -R"span[hosts=1]"'
        else:
            span = ""
        id, sources, species = ids[i], source[i], spec[i]
        if bsub_command == "Y":
            bsub = (
                "bsub -P team205 -G teichlab "
                + queue
                + " -o logs/out_SPLITTING_"
                + id
                + " -J "
                + id
                + " "
                + mem
                + " "
                + span
                + " "
            )
        if "1" in command:
            command1 = (
                "python "
                + "isotyper/Read_processing_and_quality.py "
                + dir
                + " "
                + id
                + " "
                + sample
                + " "
                + gene_types
                + " "
                + pair
                + " "
                + species
                + " "
                + sources
                + " "
                + str(100)
                + " "
                + primers
                + " "
                + platforms
                + " 1 "
                + other
                + " "
                + reverse_primer_group[i]
            )
            commands.append(bsub + command1)
        if "2" in command:
            command1 = (
                "python "
                + "isotyper/Read_processing_and_quality.py "
                + dir
                + " "
                + id
                + " "
                + sample
                + " "
                + gene_types
                + " "
                + pair
                + " "
                + species
                + " "
                + sources
                + " "
                + str(100)
                + " "
                + primers
                + " "
                + platforms
                + " 2 "
                + other
                + " "
                + reverse_primer_group[i]
            )
            commands.append(bsub + command1)
        if "3" in command:
            command1 = (
                "python "
                + "isotyper/Read_processing_and_quality.py "
                + dir
                + " "
                + id
                + " "
                + sample
                + " "
                + gene_types
                + " "
                + pair
                + " "
                + species
                + " "
                + sources
                + " "
                + str(100)
                + " "
                + primers
                + " "
                + platforms
                + " 3 "
                + other
                + " "
                + reverse_primer_group[i]
            )
            commands.append(bsub + command1)
        if "4" in command:
            command2 = (
                "python "
                + "isotyper/Generate_repertoire_statistics.py "
                + dir
                + "ORIENTATED_SEQUENCES/ANNOTATIONS/ "
                + id
                + " "
                + dir
                + "ORIENTATED_SEQUENCES/NETWORKS/Fully_reduced_"
                + id
                + ".fasta "
                + dir
                + "ORIENTATED_SEQUENCES/Filtered_ORFs_sequences_all_"
                + id
                + ".fasta "
                + gene_types
                + " "
                + species
                + " "
                + dir
                + "ORIENTATED_SEQUENCES/NETWORKS/Cluster_identities_"
                + id
                + ".txt STATISTICS "
                + reverse_primer_group[i]
            )
            commands.append(bsub + command2)
        if "ISO1" in command:
            command1 = (
                "python "
                + "isotyper/IsoTyper_1.0.py "
                + id
                + " "
                + id
                + " "
                + dir
                + " "
                + species
                + " "
                + reverse_primer_group[i]
            )
            commands.append(bsub + command1)
    # desc = file.replace("Samples_", "").replace(".txt", "")
    # if "11" in command:
    #     command11 = (
    #         "R CMD BATCH "
    #         + dirs[0]
    #         + "ORIENTATED_SEQUENCES/Network_generation_"
    #         + desc
    #         + ".R"
    #     )
    #     commands.append(bsub + command11)
    for comm in commands:
        if print_command == "Y":
            print(comm, "\n")
        if run_command == "Y":
            os.system(comm)
