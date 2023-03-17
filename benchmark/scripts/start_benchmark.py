import os

def read_constants(filename):
    constants_with_costs = []
    with open(filename, "r") as f:
        for line in f.readlines():
            line = line.replace("\n", "")
            line = line.replace("\r", "")
            if len(line) == 0:
                continue
            elements = line.split(";")
            if len(elements) < 2:
                continue
            constants_with_costs.append([elements[0], elements[1]])
    return constants_with_costs

def get_setup(bench_type):
    if bench_type == "scm_reduced":
        filename = "benchmark/scm/scm_cost_20_reduced.csv"
        name_tag = "constant"
    elif bench_type == "scm":
        filename = "benchmark/scm/scm_cost_20.csv"
        name_tag = "constant"
    elif bench_type == "mcm":
        filename = "benchmark/mcm/benchmarkfilter2d_cost.csv"
        name_tag = "counting"
    elif bench_type == "unsigned_mcm":
        filename = "benchmark/mcm/benchmarkfilter2d_cost.csv"
        name_tag = "counting"
    elif bench_type == "enumerate_unsigned":
        filename = "benchmark/mcm/benchmarkfilter2d_cost.csv"
        name_tag = "counting"
    elif bench_type == "enumerate_signed_negative":
        filename = "benchmark/mcm/benchmarkfilter2d_cost.csv"
        name_tag = "counting"
    elif bench_type == "enumerate_signed_all":
        filename = "benchmark/mcm/benchmarkfilter2d_cost.csv"
        name_tag = "counting"
    else:
        raise Exception(f"invalid benchmark specified: {bench_type}")
    return filename, name_tag

def main():
    bench_type_tuples = [("scm_lagoon_metodi",0,0), ("scm_reduced",0,0), ("scm_reduced",1,1), ("scm",0,0), ("scm",0,1), ("mcm",0,0), ("mcm",1,0), ("unsigned_mcm",1,0), ("enumerate_unsigned",0,1), ("enumerate_signed_negative",0,1), ("enumerate_signed_all",0,1)]
    for bench_type_tuple in bench_type_tuples:
        do_it(bench_type_tuple)
        bit_str = "with bit-level cost minimization" if bench_type_tuple[1] == 1 else "without bit-level cost minimization"
        shift_str = "including post-add right shifter" if bench_type_tuple[2] == 1 else "without post-add right shifter"
        print(f"finished benchmark {bench_type_tuple[0]} :)")

def do_it(bench_type_tuple):
    # setup
    bench_type = bench_type_tuple[0]
    also_minimize_full_adders = bench_type_tuple[1]
    allow_post_add_right_shift = bench_type_tuple[2]
    filename, name_tag = get_setup(bench_type)
    # paths
    result_dir = "benchmark/results"
    binary = "satscm"
    solver = "CaDiCaL"
    # timeout
    timeout_mul = 1 # standard timeout: 1h
    if "mcm" in bench_type:
        timeout_mul = 2 # 2h for mcm
    elif "lagoon" in bench_type:
        timeout_mul = 24*7 # 1 week for very hard scm instances
    elif "enumerate" in bench_type:
        timeout_mul = 24*7 # 1 week for enumeration
    timeout = 3600*timeout_mul
    if "enumerate" in bench_type:
        enumerate_all = 1
    else:
        enumerate_all = 0
    threads = 1
    quiet = 1
    write_cnf = 0
    min_num_add = 0
    allow_negative_numbers = 0 if bench_type == "unsigned_mcm" else also_minimize_full_adders
    allow_sign_inversion = 0 if bench_type == "unsigned_mcm" or allow_negative_numbers == 0 else -1
    if "enumerate" in bench_type:
        if bench_type == "enumerate_unsigned":
            allow_sign_inversion = 0
            allow_negative_numbers = 0
        elif bench_type == "enumerate_signed_negative":
            allow_sign_inversion = -1
            allow_negative_numbers = 1
        else: # bench_type == "enumerate_signed_negative"
            enumerate_signed_all = 1
            allow_negative_numbers = 1
    # create directories if they don't already exist
    if not os.path.exists(f"{result_dir}/"):
        os.mkdir(f"{result_dir}/")
    if not os.path.exists(f"{result_dir}/{bench_type}/"):
        os.mkdir(f"{result_dir}/{bench_type}/")
    if not os.path.exists(f"{result_dir}/{bench_type}/shift_{allow_post_add_right_shift}_FA_{also_minimize_full_adders}/"):
        os.mkdir(f"{result_dir}/{bench_type}/shift_{allow_post_add_right_shift}_FA_{also_minimize_full_adders}/")
    # read constants incl. costs
    constants_with_costs = read_constants(filename)
    i = 0
    for [constant, costs] in constants_with_costs:
        # check if experiment results already exist
        if name_tag == "constant":
            filename_temp = f"{constant}"
        elif name_tag == "counting":
            filename_temp = f"{i}"
            i += 1
        result_filename = f"{result_dir}/{bench_type}/shift_{allow_post_add_right_shift}_FA_{also_minimize_full_adders}/{filename_temp}.txt"
        if os.path.isfile(result_filename):
            # result already exists
            continue
        # create file
        command = f"touch {result_filename}"
        os.system(command)
        # run experiment
        fa_min_str = "without FA minimization"
        if also_minimize_full_adders:
            fa_min_str = "with FA minimization"
        print(f"C = '{constant}' expected #adders = '{costs}' {fa_min_str}")
        command = f"{binary} {constant} {solver} {timeout} {threads} {quiet} {also_minimize_full_adders} {allow_post_add_right_shift} {allow_negative_numbers} {write_cnf} {allow_sign_inversion} {min_num_add} {enumerate_all} 1>log.txt 2>>{result_filename}"
        os.system(command)
    
if __name__=="__main__":
    main()
