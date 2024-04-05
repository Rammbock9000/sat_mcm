from start_benchmark import do_it, get_num_threads_from_user

def main():
    num_threads = get_num_threads_from_user()
    worker_threads = [None for _ in range(num_threads)]
    # tuple: (<bench-type>, <subdir-name>, <min-full-adders>, <allow-right-shift>, <pipelining>)
    # if <pipelining> is not specified (i.e. len(tuple) == 4), it is assumed to be 0
    bench_type_tuples = [("pmcm","shift_1_FA_1",1,1,1)]
    for bench_type_tuple in bench_type_tuples:
        do_it(bench_type_tuple, "kissat", worker_threads)
        print(f"finished benchmark {bench_type_tuple[0]} :)")

if __name__=='__main__':
    main()
