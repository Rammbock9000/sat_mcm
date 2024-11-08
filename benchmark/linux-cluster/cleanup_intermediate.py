import os
import sys


def get_int_from_user_arg(idx, default_value):
    argc = len(sys.argv)
    if len(sys.argv) < idx+1 or len(sys.argv) < argc+1:
        return default_value
    try:
        return int(sys.argv[idx])
    except ValueError as e:
        print(f"unable to convert '{sys.argv[idx]}' to int: {repr(e)}")
    return default_value


def get_subdir_list(experiment_type, Ws=None, Ms=None, Ns=None):
    subdir_list = []
    if experiment_type == "cmm":
        if Ws is None:
            Ws = [4,6,8]
        if Ms is None:
            Ms = [2,3,4,5,6,7,8,9]
        if Ns is None:
            Ns = [2,3,4,5,6,7,8,9]
        for W in Ws:
            for M in Ms:
                for N in Ns:
                    subdir_list.append(f"{M}x{N}_{W}_bit")
    if experiment_type == "pcmm":
        if Ws is None:
            Ws = [4,6,8]
        if Ms is None:
            Ms = [2,3,4,5,6,7,8,9]
        if Ns is None:
            Ns = [2,3,4,5,6,7,8,9]
        for W in Ws:
            for M in Ms:
                for N in Ns:
                    subdir_list.append(f"{M}x{N}_{W}_bit_pipeline")
    if experiment_type == "complex":
        if Ws is None:
            Ws = [2,3,4,5,6,7,8,9]
        for W in Ws:
            subdir_list.append(f"complex_mult_{W}_bit")
    if experiment_type == "pcomplex":
        if Ws is None:
            Ws = [2,3,4,5,6,7,8,9]
        for W in Ws:
            subdir_list.append(f"complex_mult_{W}_bit_pipeline")
    if experiment_type == "conv":
        if Ws is None:
            Ws = [2,4,6,8]
        if Ns is None:
            Ns = [2,3,4,5,6,7,8,9]
        for W in Ws:
            for N in Ns:
                subdir_list.append(f"1x{N}_{W}_bit")
    if experiment_type == "pconv":
        if Ws is None:
            Ws = [2,4,6,8]
        if Ns is None:
            Ns = [2,3,4,5,6,7,8,9]
        for W in Ws:
            for N in Ns:
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
    Ws = get_int_from_user_arg(3, None)
    Ms = get_int_from_user_arg(4, None)
    Ns = get_int_from_user_arg(5, None)
    if Ws is not None:
        Ws = [Ws]
    if Ms is not None:
        Ms = [Ms]
    if Ns is not None:
        Ns = [Ns]
    subdir_list = get_subdir_list(argv[2], Ws, Ms, Ns)
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
