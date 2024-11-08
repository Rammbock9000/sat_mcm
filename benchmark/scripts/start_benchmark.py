import os
import threading
import sys
from time import sleep

def get_num_threads_from_user():
    return get_int_from_user_arg(1, 1)

def get_int_from_user_arg(idx, default_value):
    argc = len(sys.argv)
    if len(sys.argv) < idx+1 or len(sys.argv) < argc+1:
        return default_value
    try:
        return int(sys.argv[idx])
    except ValueError as e:
        print(f"unable to convert '{sys.argv[idx]}' to int: {repr(e)}")
    return default_value


class worker_thread(threading.Thread):
    def __init__(self, id, command):
        threading.Thread.__init__(self)
        self.id = id
        self.command = command

    def run(self):
        print(f"  thread id {self.id} -> executing '{self.command}'")
        os.system(self.command)

def read_constants(filename):
    constants_with_costs = []
    with open(filename, "r") as f:
        for line in f.readlines():
            line = line.replace("\n", "")
            line = line.replace("\r", "")
            if len(line) == 0:
                continue
            elements = line.split(",")
            if len(elements) < 2:
                continue
            constants_with_costs.append([elements[0], elements[1]])
    return constants_with_costs

def get_setup(bench_type, subdir=None):
    bench_type = bench_type.lower()
    if bench_type == "scm_reduced" or bench_type == "pscm_reduced":
        filename = "benchmark/inputs/scm/scm_reduced.csv"
        name_tag = "constant"
    elif bench_type == "scm":
        filename = "benchmark/inputs/scm/scm.csv"
        name_tag = "constant"
    elif bench_type == "pscm":
        filename = "benchmark/inputs/scm/scm_20.csv"
        name_tag = "constant"
    elif bench_type == "scm_lagoon_metodi":
        filename = "benchmark/inputs/scm/scm_lagoon_metodi.csv"
        name_tag = "constant"
    elif bench_type == "mcm" or bench_type == "pmcm":
        filename = "benchmark/inputs/mcm/mcm.csv"
        name_tag = "counting"
    elif bench_type == "unsigned_mcm":
        filename = "benchmark/inputs/mcm/mcm.csv"
        name_tag = "counting"
    elif bench_type == "enumerate_unsigned":
        filename = "benchmark/inputs/mcm/mcm.csv"
        name_tag = "counting"
    elif bench_type == "enumerate_signed_negative":
        filename = "benchmark/inputs/mcm/mcm.csv"
        name_tag = "counting"
    elif bench_type == "enumerate_signed_all":
        filename = "benchmark/inputs/mcm/mcm.csv"
        name_tag = "counting"
    elif bench_type == "complex_mult":
        filename = "benchmark/inputs/cmm/complex_8.csv"
        name_tag = "constant"
    elif bench_type == "fft_rotators":
        filename = "benchmark/inputs/cmm/fft_rotators.csv"
        name_tag = "constant"
    elif bench_type == "low_complexity_rotators":
        filename = "benchmark/inputs/cmm/low_complexity_rotators.csv"
        name_tag = "constant"
    elif bench_type == "fft_rot_kaya_garrido_2024":
        filename = "benchmark/inputs/cmm/fft_rot_kaya_garrido_2024.csv"
        name_tag = "constant"
    elif bench_type == "fft_rot_garrido_malagon_2021":
        filename = "benchmark/inputs/cmm/fft_rot_garrido_malagon_2021.csv"
        name_tag = "constant"
    elif bench_type == "rpag_cmm":
        filename = "benchmark/inputs/cmm/rpag_cmm.csv"
        name_tag = "constant"
    elif bench_type == "cmm_rnd":
        filename = f"benchmark/inputs/cmm_rnd/{subdir.replace('_pipeline', '')}.csv"
        name_tag = "constant"
    elif bench_type == "sop_symmetry":
        filename = "benchmark/inputs/sop/sop_symmetry_11.csv"
        name_tag = "constant"
    else:
        raise Exception(f"invalid benchmark specified: {bench_type}")
    return filename, name_tag

def main():
    bench_type_tuples = [("mcm","shift_0_FA_0",0,0), ("mcm","shift_0_FA_1",1,0), ("unsigned_mcm","shift_0_FA_1",1,0), ("scm_lagoon_metodi","shift_0_FA_0",0,0), ("scm_reduced","shift_0_FA_0",0,0), ("scm_reduced","shift_1_FA_1",1,1), ("scm","shift_0_FA_0",0,0), ("scm","shift_1_FA_0",0,1), ("enumerate_unsigned","shift_1_FA_0",0,1), ("enumerate_signed_negative","shift_1_FA_0",0,1), ("enumerate_signed_all","shift_1_FA_0",0,1)]
    num_threads = get_num_threads_from_user()
    worker_threads = [None for _ in range(num_threads)]
    for bench_type_tuple in bench_type_tuples:
        worker_threads = do_it(bench_type_tuple, "cadical", worker_threads)
        bit_str = "with bit-level cost minimization" if bench_type_tuple[1] == 1 else "without bit-level cost minimization"
        shift_str = "including post-add right shifter" if bench_type_tuple[2] == 1 else "without post-add right shifter"
    # now wait for all remaining threads to finish
    for worker in worker_threads:
        if worker is None:
            continue
        worker.join()
    print(f"finished :)")

def do_it(bench_type_tuple, solver="CaDiCaL", worker_threads=[None]):
    # setup
    min_adder_depth = 0
    equalize_output_stages = 0
    bench_type = bench_type_tuple[0]
    subdir_name = bench_type_tuple[1]
    also_minimize_full_adders = bench_type_tuple[2]
    allow_post_add_right_shift = bench_type_tuple[3]
    pipelining = 0
    if len(bench_type_tuple) > 4:
        pipelining = bench_type_tuple[4]
    if pipelining:
        equalize_output_stages = 1
        pipe_str = "with pipelining"
    else:
        pipe_str = "without pipelining"
    print(f"executing experiment '{bench_type}/{subdir_name}' ({also_minimize_full_adders}/{allow_post_add_right_shift}) {pipe_str} now")
    filename, name_tag = get_setup(bench_type, subdir_name)
    # paths
    result_dir = "benchmark/results"
    binary = "./satmcm"
    # timeout
    timeout_mul = 1 # standard timeout: 1h
    if "mcm" in bench_type:
        timeout_mul =  2   # 2h for mcm
    elif "lagoon" in bench_type:
        timeout_mul = 24*7 # 1 week for very hard scm instances
    elif "enumerate" in bench_type:
        timeout_mul = 24*7 # 1 week for enumeration
    elif "complex_mult" in bench_type or "sop_symmetry" in bench_type:
        timeout_mul = 24*5 # 5 days for complex multiplications and sop instances
    elif "fft" in bench_type or "rotators" in bench_type or "rpag_cmm" in bench_type:
        timeout_mul = 24*7 # 7 days for cmm optimization (it's hard)
    elif "cmm" in bench_type:
        timeout_mul = 12*1 # 1/2 day for "normal" cmm experiments based on random matrices
    #if pipelining:
    #    timeout_mul *= 2   # pipelining is even harder
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
            allow_sign_inversion = 1
            allow_negative_numbers = 1
    if bench_type == "complex_mult" or bench_type == "sop_symmetry":
        allow_sign_inversion = 0
        allow_negative_numbers = 1
    elif "rotators" in bench_type or "cmm" in bench_type:
        allow_sign_inversion = 1
        allow_negative_numbers = 1

    # create directories if they don't already exist
    if not os.path.exists(f"{result_dir}/"):
        os.mkdir(f"{result_dir}/")
    if not os.path.exists(f"{result_dir}/{bench_type}/"):
        os.mkdir(f"{result_dir}/{bench_type}/")
    if not os.path.exists(f"{result_dir}/{bench_type}/{subdir_name}/"):
        os.mkdir(f"{result_dir}/{bench_type}/{subdir_name}/")
    # init thread management system (a simple list lol)
    #worker_threads = [None for _ in range(num_worker_threads)]
    # read constants incl. costs
    constants_with_costs = read_constants(filename)
    i = 0
    for [constant, costs] in constants_with_costs:
        # check if experiment results already exist
        if name_tag == "constant":
            filename_temp = f"{constant}"
        else: # name_tag == "counting":
            filename_temp = f"{i}"
            i += 1
        result_filename = f"{result_dir}/{bench_type}/{subdir_name}/{filename_temp}.txt".replace(";","_")
        if os.path.isfile(result_filename):
            # result already exists
            continue
        # create file
        command = f'touch {result_filename}'
        os.system(command)
        # run experiment
        fa_min_str = "without FA minimization"
        if also_minimize_full_adders:
            fa_min_str = "with FA minimization"
        # wait until we can spawn a new worker thread
        worker_id = -1
        while worker_id < 0:
            for potential_worker_id, worker in enumerate(worker_threads):
                if worker is None:
                    # found uninitialized worker
                    worker_id = potential_worker_id
                    break
                if worker.is_alive():
                    # found running worker
                    continue
                # found idle (finished) worker
                worker.join()
                worker_id = potential_worker_id
            sleep(0.1)
        # start new worker
        command = f'{binary} "{constant}" solver_name={solver} timeout={timeout} threads={threads} quiet={quiet} minimize_full_adders={also_minimize_full_adders} allow_post_adder_right_shift={allow_post_add_right_shift} allow_negative_coefficients={allow_negative_numbers} write_cnf_files={write_cnf} allow_coefficient_sign_inversion={allow_sign_inversion} min_num_adders={min_num_add} enumerate_all={enumerate_all} min_adder_depth={min_adder_depth} pipelining={pipelining} equalize_output_stages={equalize_output_stages} 1>log_{worker_id}.txt 2>>{result_filename}'
        worker_threads[worker_id] = worker_thread(worker_id, command)
        worker_threads[worker_id].start()
        sleep(0.1)
    return worker_threads

    
if __name__=="__main__":
    main()
