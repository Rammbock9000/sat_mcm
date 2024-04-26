import os
import sys


def get_subdir_list(experiment_type):
    subdir_list = []
    if experiment_type == "cmm":
        for W in [4,6,8]:
            for M in [2,3,4,5,6,7,8,9]:
                for N in [2,3,4,5,6,7,8,9]:
                    subdir_list.append(f"{M}x{N}_{W}_bit")
    if experiment_type == "pcmm":
        for W in [4,6,8]:
            for M in [2,3,4,5,6,7,8,9]:
                for N in [2,3,4,5,6,7,8,9]:
                    subdir_list.append(f"{M}x{N}_{W}_bit_pipeline")
    if experiment_type == "complex":
        for W in [2,3,4,5,6,7,8,9]:
            subdir_list.append(f"complex_mult_{W}_bit")
    if experiment_type == "pcomplex":
        for W in [2,3,4,5,6,7,8,9]:
            subdir_list.append(f"complex_mult_{W}_bit_pipeline")
    if experiment_type == "conv":
        for W in [2,4,6,8]:
            for N in [2,3,4,5,6,7,8,9]:
                subdir_list.append(f"1x{N}_{W}_bit")
    if experiment_type == "pconv":
        for W in [2,4,6,8]:
            for N in [2,3,4,5,6,7,8,9]:
                subdir_list.append(f"1x{N}_{W}_bit_pipeline")
    return subdir_list


def main():
    argv = sys.argv
    argc = len(argv)
    if argc < 2:
        raise Exception("need program arguments: path to directory for cleanup AND experiment type")
    if argc < 3:
        raise Exception("need program argument: experiment type")
    basedir = argv[1]
    subdir_list = get_subdir_list(argv[2])
    for subdir in subdir_list:
        dir_path = f"{basedir}/{subdir}"
        if not os.path.isdir(dir_path):
            print(f"skip missing directory '{dir_path}'")
            continue
        for file in os.listdir(os.fsencode(dir_path)):
            filename = os.fsdecode(file)
            if filename.startswith("._"):
                continue
            if not filename.endswith(".txt"):
                continue
            filename = f"{dir_path}/{filename}"
            # check if experiment ran to completion
            found_add_opt_line = False
            found_fas_opt_line = False
            with open(filename, "r") as f:
                for line in f:
                    line = line.rstrip()
                    if "#Add optimal = " in line:
                        found_add_opt_line = True
                    if "#FAs optimal = " in line:
                        found_fas_opt_line = True
            if found_add_opt_line and found_fas_opt_line:
                continue
            # delete file if not
            print(f"deleting file '{filename}'")
            os.remove(filename)


if __name__ == '__main__':
    main()
