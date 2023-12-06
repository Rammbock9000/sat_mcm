from start_benchmark import do_it, get_num_threads_from_user

def main():
    num_threads = get_num_threads_from_user()
    # tuple: (<bench-type>, <subdir-name>, <min-full-adders>, <allow-right-shift>, <pipelining>)
    # if <pipelining> is not specified (i.e. len(tuple) == 4), it is assumed to be 0
    bench_type_tuples = [("sop_symmetry","shift_1_FA_0",0,1)]
    for bench_type_tuple in bench_type_tuples:
        do_it(bench_type_tuple, "kissat", num_worker_threads=num_threads)
        print(f"finished benchmark {bench_type_tuple[0]} :)")

if __name__=='__main__':
    main()