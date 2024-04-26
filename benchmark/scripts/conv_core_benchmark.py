from start_benchmark import do_it, get_num_threads_from_user

def main():
    num_threads = get_num_threads_from_user()
    worker_threads = [None for _ in range(num_threads)]
    # tuple: (<bench-type>, <subdir-name>, <min-full-adders>, <allow-right-shift>, <pipelining>)
    # if <pipelining> is not specified (i.e. len(tuple) == 4), it is assumed to be 0
    #template: bench_type_tuples = [("cmm_rnd","subdir",0,1,0)]
    bench_type_tuples = []
    # add experiments "1xN vs W bits"
    for W in [2,4,6,8]:
        for N in [2,3,4,5,6,7,8,9]:
            bench_type_tuples.append(("cmm_rnd", f"1x{N}_{W}_bit",1,1,0))
    # run benchmarks
    for bench_type_tuple in bench_type_tuples:
        worker_threads = do_it(bench_type_tuple, "kissat", worker_threads)
    # now wait for all remaining threads to finish
    for worker in worker_threads:
        if worker is None:
            continue
        worker.join()
    print(f"finished :)")


if __name__=='__main__':
    main()
