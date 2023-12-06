from start_benchmark import do_it

def main():
    # tuple: (<bench-type>, <subdir-name>, <min-full-adders>, <allow-right-shift>, <pipelining>)
    # if <pipelining> is not specified (i.e. len(tuple) == 4), it is assumed to be 0
    bench_type_tuples = [("complex_mult","shift_1_FA_0_inv_1",0,1)]
    for bench_type_tuple in bench_type_tuples:
        do_it(bench_type_tuple, "kissat", num_worker_threads=1)
        print(f"finished benchmark {bench_type_tuple[0]} :)")

if __name__=='__main__':
    main()