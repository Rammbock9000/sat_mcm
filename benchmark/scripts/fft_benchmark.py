from start_benchmark import do_it, get_num_threads_from_user, get_int_from_user_arg


def main():
    num_threads = get_num_threads_from_user()
    worker_threads = [None for _ in range(num_threads)]
    # tuple: (<bench-type>, <subdir-name>, <min-full-adders>, <allow-right-shift>, <pipelining>)
    # if <pipelining> is not specified (i.e. len(tuple) == 4), it is assumed to be 0
    #template: bench_type_tuples = [("cmm_rnd","subdir",0,1,0)]
    worker_threads = do_it(("fft_rot_garrido_malagon_2021", f"fft_rot_garrido_malagon_2021-pipe_0",1,1,0), "kissat", worker_threads)
    worker_threads = do_it(("fft_rot_garrido_malagon_2021", f"fft_rot_garrido_malagon_2021-pipe_1",1,1,1), "kissat", worker_threads)
    # now wait for all remaining threads to finish
    for worker in worker_threads:
        if worker is None:
            continue
        worker.join()
    print(f"finished :)")


if __name__=='__main__':
    main()
